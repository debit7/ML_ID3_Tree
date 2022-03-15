import pandas as pd
import math
from anytree import Node, RenderTree

with open("fishing.data", "r") as input:
    whole_data = input.read().split("\n\n") 
Heading=whole_data[0].replace('#','')
def Get_Attributes(Adata,Cdata):
    x=Adata.split('\n')
    attributes=[]
    for y in x:
        atts=y.split('#')[1].split(':')[0]
        attributes.append(atts)
    attributes.pop(0)
    c=Cdata.split('\n')
    cattributes=[]
    for y in c:
        atts=y.split('#')[1].split(':')[0]
        cattributes.append(atts)
    cattributes.pop(0)
    attributes.append(cattributes[0])
    return attributes
def Get_Data(Rdata):
    Rdata=Rdata.split('\n')
    Rdata.pop(0)
    lst = []
    for ele in Rdata:
        line = ele.replace(' ','').split(',')
            
        lst.append(line)
    # lst.pop(-1)
    return lst
def Get_Targets(Tdata):
    Tdata=Tdata.split(':')
    targets=Tdata[1].split(',')

    return targets
def Get_Attribute_classes(Cdata):
    Data=Cdata.split('\n')
    Data.pop(0)
    set_data=[]
    for dta in Data:
        new=dta.split(':')
        gg=new[1].split(',')
        set_data.append(gg)
    return set_data


Targets=Get_Targets(whole_data[2])
# print(Targets)
Headers=Get_Attributes(whole_data[1],whole_data[2])
print(Headers)
df = pd.DataFrame(Get_Data(whole_data[3]),columns =Headers) 
df[df.columns] = df.apply(lambda x: x.str.strip())
Atts_labels=Get_Attribute_classes(whole_data[1])
print(Atts_labels)
heads=Get_Attributes(whole_data[1],whole_data[2])
heads.pop(len(heads)-1)
Headers_labels_mapping= {'Attributes':heads,'Labels':Atts_labels}
Headers_labels_mapping=pd.DataFrame(Headers_labels_mapping)
print(Headers_labels_mapping)
def count(data,colname,label,target):
    condition = (data[colname] == label) & (data[Headers[-1]] == target)
    return len(data[condition])
def Get_Probability_Total(Df,Headers,Targets):
    total_probs=[]
    total_counts=[]
    
    for trt in Targets:

        
        total_counts.append(count(Df,Headers[-1],trt.strip(),trt.strip()))
    total_outcomes=sum(total_counts)
    
    for prob in total_counts:
        total_probs.append(prob/total_outcomes)
    
    # target_probs = dict(zip(Targets, total_probs))
    return total_probs
print(Get_Probability_Total(df,Headers,Targets))
# print (Atts_labels)
def Calculate_entropy(List_of_probs):
    

    Entropy=0
    for probs in List_of_probs:
        if probs==0:
            return 0
        Entropy+=probs*math.log2(probs)*(-1)
    return Entropy
Total_entropy=Calculate_entropy(Get_Probability_Total(df,Headers,Targets))
def Calculate_Gain(Total_entropy,Ratio,Entropy):
    i=0
    Gain=0
    
    for ent in Entropy:
        Gain+=Ratio[i]*ent
        i+=1
    Final_Gain=Total_entropy-Gain
    return Final_Gain
    



def Countings(Atts_labels):
    countings=[]
    gain_df=[]
    i=0
    for labels in Atts_labels:
            entropies=[]
            ratios=[]
            for dta in labels:
                j=0
                
                prob_arrays=[]
                totals=0
                s=0   
                while(s<len(Targets)):
                        totals+=count(df,Headers[i],dta.strip(),Targets[s].strip())
                        prob_arrays.append(count(df,Headers[i],dta.strip(),Targets[s].strip()))
                        s=s+1
                entropy_ars=[]
                for ars in prob_arrays:
                    entropy_ars.append(ars/totals)
                real_entropy=Calculate_entropy(entropy_ars)
                entropies.append(real_entropy)
                ratios.append(sum(prob_arrays)/len(df))
                countings.append([Headers[i],dta.strip(),prob_arrays,real_entropy])
            Gain=Calculate_Gain(Total_entropy,ratios,entropies)
            gain_df.append([Headers[i],Gain])


            
                
                # while(j<len(Targets)):
                #     print(Headers[i],dta.strip(),Targets[j])
                    
                    
                        
                #     countings.append([Headers[i],dta.strip(),Targets[j].strip(),count(df,Headers[i],dta.strip(),Targets[j].strip()),totals,prob_arrays,real_entropy])
                #     j+=1
            i+=1
    # print(countings)
    return pd.DataFrame(countings,columns=['Attributes','Labels','Array','entropy']),pd.DataFrame(gain_df,columns=['Attributes','Gain'])

Calculation_Frame,Gain_dataframe=Countings(Atts_labels)
# Calculation_Frame.groupby(['Attributes','Labels']).apply(lambda x: Calculate_entropy(x['entropy']))
# Calculation_Frame.groupby(['Attributes','Labels']).apply(lambda x: Calculate_entropy(x['value1'])).astype(int)
        
print(Gain_dataframe)

def Find_Max(df):
    col = "Gain"
    max_x = df.loc[df[col].idxmax()]
    return max_x['Attributes']

def Draw_Tree_Diagram(Headers_labels_mapping,Highest_Gain):
    # Headers_labels_mapping.where(Headers_labels_mapping['Attributes']==Highest_Gain, inplace=True)
    rslt_df = Headers_labels_mapping[Headers_labels_mapping['Attributes']==Highest_Gain]
    labels=rslt_df['Labels']
    udo = Node(Highest_Gain)
    for lst in labels:
        for ss in lst:
            ss = Node(ss, parent=udo)
    for pre, fill, node in RenderTree(udo):
        print("%s%s" % (pre, node.name))

    
    
    

Draw_Tree_Diagram(Headers_labels_mapping,Find_Max(Gain_dataframe))
