from itertools import permutations
import glob,argparse
from basic_modules import *
import pandas as pd
import matplotlib.pyplot as plt

### read args ###

parser = argparse.ArgumentParser(description='input for running go enrichment analysis')
parser.add_argument('time1',type=str, help='the first timestamp, in the format of %Y-%m-%d')
parser.add_argument('time2',type=str, help='the second timestamp, in the format of %Y-%m-%d')
parser.add_argument('data_path',type=str, help='path to the data')
parser.add_argument('--visualize', help='if the enrichment comparsion will be visualized')

args=parser.parse_args()

data_path=args.data_path   #'./example_data/covidGroups/'
time1=args.time1 #'2017-09-10'
time2=args.time2 #'2023-09-10'
visualization=args.visualize
groups={}

def main(data_path=data_path,time1=time1,time2=time2,visualization=visualization):
    for path in glob.glob(data_path+'*'): ## retrieve all groups 
        group_number=path.split('group')[1]
        groupname='group'+str(group_number)
        groupname=groupname.split('.')[0]
        genes=[]
        with open(path,'r') as f:
            for line in f.readlines():
                genes.append(line.rstrip())
        groups[groupname]=genes

    f1=goenrichment(time1,anno_count=True,ontology_desc=True) ## initialize
    f2=goenrichment(time2,anno_count=True,ontology_desc=True) ## initialize

    records=pd.DataFrame(columns=['group','Jac_index'])
    n=0

    for x,y in groups.items():
        result1=f1.enrichment(y).sort_values(by='scores_adj').head(20)
        result2=f2.enrichment(y).sort_values(by='scores_adj').head(20)
        result1.to_csv('./output/data/{0}_{1}.csv'.format(time1.split('-')[0]+time1.split('-')[1],x))
        result2.to_csv('./output/data/{0}_{1}.csv'.format(time2.split('-')[0]+time2.split('-')[1],x))
        head1=result1.GOterms.to_list()
        head2=result2.GOterms.to_list()
        share=set(head1)&set(head2)
        h1_only=set(head1)-set(head2)
        h2_only=set(head2)-set(head1)
        total=share|h1_only|h2_only
        records.loc[n]=[x,len(share)/len(total)]
        n+=1
        if visualization:
            if len(share)<20:
                pairs1=permutations(head1,2)
                pairs2=permutations(head2,2)
                edges1=[]
                edges2=[]
                for pair in pairs1:
                    try:
                        sp=nx.shortest_path(f1.ontology,pair[0],pair[1])
                    except:
                        continue
                    else:
                        if len(sp)<=4:
                            path_edges = [(sp[n], sp[n+1]) for n in range(len(sp)-1)]
                            edges1.extend(path_edges)
                for pair in pairs2:
                    try:
                        sp2=nx.shortest_path(f2.ontology,pair[0],pair[1])
                    except:
                        continue
                    else:
                        if len(sp2)<=4:
                            path_edges2 = [(sp2[n], sp2[n+1]) for n in range(len(sp2)-1)]
                            edges2.extend(path_edges2)
                edge_list=set(edges2)|set(edges1)
                G=nx.Graph()
                G.add_nodes_from(total)
                G.add_edges_from(edge_list)
                nodes=G.nodes
                color = ['grey' if x in share else 'red' if x in h1_only else 'blue' if x in h2_only else 'white' for x in nodes]
                plt.figure(figsize=(8, 6))
                pos=nx.spring_layout(G)
                nx.draw_networkx_nodes(G,node_color=color,pos=pos)
                nx.draw_networkx_edges(G,edge_color='grey',pos=pos)
                nx.draw_networkx_labels(G,font_size=6,pos=pos)
                plt.savefig("./output/figure/Group{0}{1}&{2}.png".format(x,time1,time2), format="PNG")
                plt.close()
        else:
            continue
    print('complete')






if __name__=='__main__':
    main(data_path=data_path,time1=time1,time2=time2,visualization=visualization)