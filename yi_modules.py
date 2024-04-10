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


r=requests.get('https://go-data-product-release.s3.amazonaws.com/?list-type=2&delimiter=/&prefix=')
pattern=re.compile('>[0-9]{4}\-[0-9]{2}\-[0-9]{2}/')
all_date=pattern.findall(r.text)
for x in all_date:
    all_date.remove(x)
    x=x.replace('>','')
    x=x.replace('/','')
    all_date.insert(0,x)
all_date.sort()



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
            annotation,ontology=corresponding_time(time,taxo)
            if re.match(r'^http:/{2}\w.+$',annotation):
                output=requests.get(annotation,stream=True)
                total=int(output.headers.get('content-length'),0)
                with open('annotation_go.gz','wb') as file,tqdm.tqdm(desc='annotation_go.gz',total=total,unit='iB',unit_scale=True,unit_divisor=1024) as bar:
                    for data in output.iter_content(chunk_size=1024):
                        size=file.write(data)
                        bar.update(size)
                f1 = gzip.open('annotation_go.gz', 'rb')

                items=[]
                name=['gene product','relationship','GO term','evidence code']
                for line in f1.readlines():
                    line=line.decode('utf-8')
                    if line.startswith('!')==False:
                        elements=line.split('\t')
                        items.append((elements[1],elements[3],elements[4],elements[6]))
                f1.close()
                os.remove('annotation_go.gz')
                annotation_f=pd.DataFrame.from_records(items,columns=name)
            
        else:
            print(args)
            annotation_local=args[0]+'goa_human.gaf'
            with open(annotation_local, 'r') as out:
                num = 0
                for line1 in out.readlines():
                    if line1.startswith('!'):
                        num += 1
                        continue
                    else:
                        break
            annotation_f = pd.read_csv(annotation_local, sep='\t', skiprows=num, header=None,usecols=[1,3,4,6],\
                                       names=['gene product','relationship','GO term','evidence code'])
        #G_a = nx.from_pandas_edgelist(annotation_f, source='gene product', target='GO term', \
        #                              edge_key='relationship', create_using=nx.MultiDiGraph())
            ontology=args[0]+'go.obo'
        G_o = obonet.read_obo(ontology)
        #G_o_obsolete=obonet_modif(ontology)
        # annotation_gene=annotation_f[['gene product','GO term']].copy()
        # annotation_gene.drop_duplicates(inplace=True)
        # G_a_gene=nx.from_pandas_edgelist(annotation_gene,source='gene product',target='GO term',\
        #                                  create_using=nx.MultiDiGraph())

        if extra_r==False:
            G_o.remove_edges_from([(x,y) for x,y,z in G_o.edges(keys=True) if (z in ['is_a','part_of'])==False])
        
        self.ontology=G_o
        self.annotation=annotation_f
        if ontology_desc:
            node_map={y:x for x,y in enumerate(G_o.nodes())} ##GOID:digits
            node_back={x:y for x,y in enumerate(G_o.nodes())} ##digits:GOID
            dig_go=nx.relabel_nodes(G_o,node_map) ##network with digits
            #inverse={node_back[x]:list(nx.ancestors(dig_go,x)) for x in node_back.keys()} ##GOID:digits
            real={x:set(nx.ancestors(dig_go,x))|{x} for x in node_back.keys()}
            #real={x:[node_back[y1] for y1 in y] for x,y in inverse.items()} ##GOID:GOID_list
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
            
        #self.genes=G_a_gene
        #self.ontology_with_obsolete=G_o_obsolete
        biological_process = self.categories()[0]
        molecular_function = self.categories()[1]
        cellular_component = self.categories()[2]
        self.aspects = ['biological_process', 'molecular_function', 'cellular_component']
        anno_nums=[]
        for aspect in self.aspects:
            asp_set=eval(aspect)
            anno_num=len(self.annotation.loc[self.annotation['GO term'].isin(asp_set)])
            anno_nums.append(anno_num)
        self.anno_nums=anno_nums
        self.n_ns={node: data['namespace'] for node, data in self.ontology.nodes(data=True)}
       
    def categories(self):
        G_o = self.ontology
        bp = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'biological_process']
        mf = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'molecular_function']
        cc = [x for x in G_o.nodes() if G_o.nodes[x]['namespace'] == 'cellular_component']
        return [bp, mf, cc]
    
    def MICA(self,pairs):
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
    
    

    def Resnik_IC(self,entities,namespace):
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

        ###1) 计算出每个subgraph的annotation数量###
        ###2) resnik 计算公式：-log(p);p=term_annotation/total_annotation
        ###3) resnik normalization -log(p)/logN<==>-log(p)/-log(1/N)
        ##step 1
        root_IC=np.log(root_num)
        ##step 2 first, find how many annotations it has

        for entity in entities_digit:
            # try:
            #     child=[subterm for subterm in onto_desc[entity]]
            # except nx.exception.NetworkXError:
            #     IC_file[entity] = anno_count(entity)#annotation.in_edges(entity))

            # else:
            #     child.append(entity)
            #     ICs = np.sum([anno_count[child_term] for child_term in child])
            ICs=self.anno_count[entity]
        ##step 2 second, compute the -log(p)
            IC=[-np.log(ICs/root_num) if ICs!=0 else 0]
    ##step 3 normalization
            IC_norm=IC/root_IC
            IC_file[entity]=IC_norm  # 这一步用于计算information content

        return IC_file

    def Seco_IC(self,entities):
        ontology=self.ontology
        aspects=['biological_process','molecular_function','cellular_component']
        maxnodes={}
        for aspect in aspects:
            maxnode=len([x for x in ontology.nodes() if ontology.nodes()[x]['namespace']==aspect])
            maxnodes[aspect]=maxnode
    
        seco_ICs={}
        for entity in entities:
            maxnode=maxnodes[ontology.nodes()[entity]['namespace']]
            hypox=len(nx.descendants(ontology,entity))
            seco_IC=1-np.log(hypox+1)/np.log(maxnode)
            seco_ICs[entity]=seco_IC
        return seco_ICs


#####load embedding data
def embedding_load(emd_path,entity_path):
    embeddings=np.load(emd_path)
    entity_index=pd.read_csv(entity_path,sep='\t',header=None)
    entity_index.columns=['index','entity']
    entity_index.set_index('index')
    index_dict=entity_index.to_dict()['entity']
    return [embeddings,index_dict]


def mint_file_load(path):
    file=pd.read_csv(path,sep='\t',header=None,usecols=[0,1,9,10,11,13,14],names=['id1','id2','taxo1','taxo2','interaction type','interaction identifier','confidence score'])
    return file 


##############interpro download#####################

def interpro_api(endpoints):
    url_base='https://www.ebi.ac.uk/interpro/api/'
    if 'cursor' in endpoints:
        url=url_base+endpoints+'&page_size=200'
    else:
        url=url_base+endpoints+'/?page_size=200'
    
    next=url
    last_page=False
    attempts=0

    output_df=pd.DataFrame(columns=['uniprotid','taxid'])
    annotation_dict={}
    n=0
    while next:
        #print(next)
        try:
            response=requests.get(next,headers={"Accept": "application/json"})
            if response.status_code!=200:
                #print(response.status_code)
                #print(next)
                break
            output=response.json()
            next=output['next']
            if not next:
                last_page=True
        except HTTPError as e:
            if e.code==408:
                sleep(10)
                continue
            else:
                if attempts<3:
                    attempts+=1
                    sleep(10)
                    continue
                else:
                    raise e
        for result in output['results']:
            uniprotID=result['metadata']['accession']
            taxid=result['metadata']['source_organism']['taxId']
            annotation_info=[f"{','.join([entry['accession'],entry['entry_type']])}" if entry is not None else '' for entry in result['entry_subset'] ]
            annotation_info=f"{';'.join(annotation_info)}"
            output_df.loc[n]=[uniprotID,taxid]
            annotation_dict[uniprotID]=annotation_info
            n+=1
    return output_df,annotation_dict,next

#########in case 503 error###############
class interpro:

    def __init__(self,endpoints=None):
        if endpoints==None:
            url='protein/reviewed/entry/pfam/taxonomy/uniprot/9606'
        else:
            url=endpoints
        self.url=url

    def interpro_api_continue(self):
        uniprot_sum=pd.DataFrame(columns=['uniprotid','taxid'])
        annotation_dict_sum={}
        url=self.url
        while len(uniprot_sum)<19000:
            test1=interpro_api(url)
            uniprot_sum=pd.concat([uniprot_sum,test1[0]],ignore_index=True)
            annotation_dict_sum.update(test1[1])
            try:
                api_1=test1[2].split('api/')[1]
                api_2=api_1.split('&page')[0]
                url=api_2
            except AttributeError:
                if test1[2] is None:
                    break
                else:
                    print(test1[2])
                    print(len(test1[0]))
                    break
        return uniprot_sum,annotation_dict_sum

    def output(self,function='family'):
        out={}
        uniprot_info,annotations=self.interpro_api_continue()
        protein_has_family=[x for x,y in annotations.items() if function in y]
        for x in protein_has_family:   
            annotation=annotations[x].split(';')
            t=[re.findall('PF[0-9]+',x)[0] for x in annotation if function in x]

            out[x]=t

        return out
    
def retrieving_embedding(embedding_vectors,terms,agg_method='mean'):
    original_vector=np.zeros(embedding_vectors.vector_size)
    for x in terms:
        vector=embedding_vectors[x]
        original_vector=np.vstack((original_vector,vector))
    if agg_method=='mean':
        original_vector=original_vector.mean(0)
    elif agg_method=='sum':
        original_vector=original_vector.sum(0)
    return original_vector



def AUC_score(label,predictions):
    fpr, tpr, thresholds = metrics.roc_curve(label, predictions)
    return metrics.auc(fpr,tpr)

class goenrichment(network):
        def __init__(self, time:datetime = None, extra_r=False, taxo='human', anno_count=False, ontology_desc=False, *args):
            super().__init__(time, extra_r, taxo, anno_count, ontology_desc, *args)
            #existing_genes=list(self.annotation['gene product'].unique())
            self.go_desc={x:set(list(nx.ancestors(self.ontology,x))+[x]) for x in self.ontology.nodes}
            annotation=self.annotation
            self.annotation_per_go={}
            for x,y in annotation.groupby('GO term')['gene product']:
                self.annotation_per_go[x]=y.to_list()
            #self.anno_per_GO={x:list(annotation.loc[annotation['GO term'].isin(go_desc[x]),'gene product'].unique()) for x in go_desc.keys()}
            self.annotation_withdesc={}
            for x,y in self.go_desc.items():
                lists = []  # 确保lists是一个已经初始化的空列表
                for y1 in y:
                    if y1 in self.annotation_per_go.keys():
                        lists.extend(self.annotation_per_go[y1])
                self.annotation_withdesc[x]=len(set(lists))
                        

        def enrichment(self,genelist,adjust_method=None):
            genes=genelist
            counts={x:set(genes)&set(y) for x,y in self.annotation_per_go.items() if len(set(genes)&set(y))>0} #GO term: gene in genelist that are annotated by this term
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