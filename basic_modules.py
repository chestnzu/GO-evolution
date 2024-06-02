import requests,re,gzip,obonet,os,re
import pandas as pd
import networkx as nx
import numpy as np
import statsmodels.stats
import statsmodels.stats.multitest
import datetime,random,tqdm
from time import sleep
from urllib.error import HTTPError
from sklearn import metrics
from scipy.stats import hypergeom
import statsmodels

def aspect(x):
    if x=='biological_process':
        return 'bp'
    elif x=='molecular_function':
        return 'mf'
    else:
        return 'cc'



###find all available version of GOKB ###
r=requests.get('https://go-data-product-release.s3.amazonaws.com/?list-type=2&delimiter=/&prefix=')
pattern=re.compile('>[0-9]{4}\-[0-9]{2}\-[0-9]{2}/')
all_date=pattern.findall(r.text)
for x in all_date:
    all_date.remove(x)
    x=x.replace('>','')
    x=x.replace('/','')
    all_date.insert(0,x)
all_date.sort()


### find the closest GO file based on the timestamp that users provided ##
def corresponding_time(time=None,taxo='human'):

    annotation_online_path='http://release.geneontology.org/0000-00-00/annotations/goa_{0}.gaf.gz'.format(taxo)
    ontology_online_path='http://release.geneontology.org/0000-00-00/ontology/go-basic.obo'
    if time==None:
        annotation_online_path='http://current.geneontology.org/annotations/goa_{0}.gaf.gz'.format(taxo)
        ontology_online_path='http://purl.obolibrary.org/obo/go-basic.obo'
        return [annotation_online_path, ontology_online_path]
    for x in all_date:
        if time<x:
            index=all_date.index(x)
            break
        else:
            index=-1
    if index==0:
        time1=all_date[0]
    else:
        time1=all_date[index-1]
    annotation_online_path = annotation_online_path.replace('0000-00-00', time1)
    ontology_online_path = ontology_online_path.replace('0000-00-00', time1)
    return [annotation_online_path, ontology_online_path]




######################construct network####################
class network:
    def __init__(self,time:datetime=None,extra_r=False,taxo='human',anno_count=False,ontology_desc=False,*args):
        if time!='none':
            annotation,ontology=corresponding_time(time,taxo) ### retrieving URL of the annotation file & ontology file based on the timepoint
            if re.match(r'^http:/{2}\w.+$',annotation):
                output=requests.get(annotation,stream=True) ### download the annotation file 
                total=int(output.headers.get('content-length'),0)
                with open('annotation_go.gz','wb') as file,tqdm.tqdm(desc='annotation_go.gz',total=total,unit='iB',unit_scale=True,unit_divisor=1024) as bar:
                    for data in output.iter_content(chunk_size=1024):
                        size=file.write(data)
                        bar.update(size)
                f1 = gzip.open('annotation_go.gz', 'rb')

                items=[]
                for line in f1.readlines():
                    line=line.decode('utf-8')
                    if line.startswith('!')==False: ## skip annotation lines 
                        elements=line.split('\t')
                        items.append((elements[1],elements[3],elements[4],elements[6]))
                f1.close()
                os.remove('annotation_go.gz')
                name=['gene product','relationship','GO term','evidence code']
                annotation_f=pd.DataFrame.from_records(items,columns=name)
            else:
                return 'invalid annotation URL'
        G_o = obonet.read_obo(ontology)

        if extra_r==False:
            G_o.remove_edges_from([(x,y) for x,y,z in G_o.edges(keys=True) if (z in ['is_a','part_of'])==False])
        
        self.ontology=G_o
        self.annotation=annotation_f
        if ontology_desc:  ## if the user want a dictionary, in which the parent term is the key, and the value is a list includes all its descendants
            node_map={y:x for x,y in enumerate(G_o.nodes())} ## GOID:digits
            node_back={x:y for x,y in enumerate(G_o.nodes())} ## digits:GOID
            dig_go=nx.relabel_nodes(G_o,node_map) ## network with digits
            real={x:set(nx.ancestors(dig_go,x))|{x} for x in node_back.keys()}
            self.node_map=node_map
            self.node_back=node_back
            self.ontology_desc=real   
        if anno_count:
            goinanno=self.annotation['GO term'].to_list()
            for x in goinanno: 
                if x not in node_map.keys():
                    node_map[x]=len(node_map)+1
                    node_back[len(node_map)]=x
            anno_count={node_map[x]:len(y) for x,y in self.annotation.groupby('GO term')['gene product']} ##this is specifically for IC calculation
            anno_count0={node_map[x]:0 for x in G_o.nodes() if node_map[x] not in anno_count.keys()}
            anno_count.update(anno_count0)
            self.anno_count={}
            for x,y in self.ontology_desc.items(): #GO term and its descendants
                count_y=0
                for item in y:
                    count_y+=anno_count[item]
                self.anno_count[x]=count_y
        biological_process = self.categories()[0]
        molecular_function = self.categories()[1]
        cellular_component = self.categories()[2]
        self.aspects = ['biological_process', 'molecular_function', 'cellular_component']
        anno_nums=[]
        for aspect in self.aspects:
            asp_set=eval(aspect)
            anno_num=len(self.annotation.loc[self.annotation['GO term'].isin(asp_set)])
            anno_nums.append(anno_num) ## the number of annotations that root terms have
        self.anno_nums=anno_nums 
        self.n_ns={node: data['namespace'] for node, data in self.ontology.nodes(data=True)} ## node and its namespace 
       
    def categories(self):
        G_o = self.ontology
        bp = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'biological_process']
        mf = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'molecular_function']
        cc = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'cellular_component']
        return [bp, mf, cc]
    
    def MICA(self,pairs): ## find the most informative common ancestors of a given pair of ontology terms
        ancestor_cache = {}
        MICA_record={}
        for v,w in pairs:
                if v not in ancestor_cache:
                    ancestor_cache[v] = nx.descendants(self.ontology, v)
                    ancestor_cache[v].add(v)
                if w not in ancestor_cache:
                    ancestor_cache[w] = nx.descendants(self.ontology, w)
                    ancestor_cache[w].add(w)

                common_ancestors = list(ancestor_cache[v] & ancestor_cache[w])

                if common_ancestors:
                    MICA_record[(v,w)]=common_ancestors
        return MICA_record
    
    

    def Resnik_IC(self,entities,namespace):  ## using Resnik method to compute the information content value of a given list of GO terms 
        n_ns=self.n_ns
        anno_nums=self.anno_nums
        aspects=self.aspects
        entities_digit=[self.node_map[x] for x in entities]
        if self.anno_count:
            anno_count=self.anno_count
        else:
            raise ValueError ('anno_count should be set as True')
        if self.ontology_desc:
            onto_desc=self.ontology_desc
        else:
            raise ValueError ('ontology_IC should be set as True')
        IC_file={}
        root_num=anno_nums[aspects.index(namespace)]
        root_IC=np.log(root_num)

        for entity in entities_digit:
            ICs=self.anno_count[entity]
            IC=[-np.log(ICs/root_num) if ICs!=0 else 0]
            IC_norm=IC/root_IC
            IC_file[entity]=IC_norm 

        return IC_file

class goenrichment(network):  ### Gene Ontology Enrichment analysis using the version of GOKB. 
        def __init__(self, time:datetime = None, extra_r=False, taxo='human', anno_count=False, ontology_desc=False, *args):
            super().__init__(time, extra_r, taxo, anno_count, ontology_desc, *args)
            self.go_desc={x:set(list(nx.ancestors(self.ontology,x))+[x]) for x in self.ontology.nodes}
            annotation=self.annotation
            self.annotation_per_go={}
            for x,y in annotation.groupby('GO term')['gene product']:
                self.annotation_per_go[x]=y.to_list()
            self.annotation_withdesc={}
            for x,y in self.go_desc.items():
                lists = [] 
                for y1 in y:
                    if y1 in self.annotation_per_go.keys():
                        lists.extend(self.annotation_per_go[y1])
                self.annotation_withdesc[x]=len(set(lists))
                        

        def enrichment(self,genelist,adjust_method=None): ## now only one adjust_method is supported 
            genes=genelist
            counts={x:set(genes)&set(y) for x,y in self.annotation_per_go.items() if len(set(genes)&set(y))>0}  #GO term: gene in genelist that are annotated by this term
            go_desc=self.go_desc
            scores=[]
            GOterm=[]
            existing_genes=list(self.annotation['gene product'].unique())
            N=len(existing_genes)
            n=len(genelist)
            Ms=[]
            ks=[]
            for x,y in go_desc.items():
                save=[]
                M=self.annotation_withdesc[x]
                Ms.append(M)
                for y1 in y: 
                    if y1 in counts.keys():
                        save.extend(list(counts[y1]))  
                GOterm.append(x)
                k=len(set(save))
                scores.append(hypergeom.sf(k-1,N,M,n))
                ks.append(k)
            scores_adj=statsmodels.stats.multitest.fdrcorrection(scores)[1]
            results=pd.DataFrame.from_dict({'GOterms':GOterm,'scores':scores,'scores_adj':scores_adj,'k':ks,'M':Ms})
            return results
