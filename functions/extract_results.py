#function to organise results into a dictionary of compartments, each compartment containing a dictionary of aggregation states, of wich contain number of particles per size fraction


def extract_by_comp(Results,compartmentNames_list):
    particle_compartmentCoding = dict(
    zip(compartmentNames_list, list(range(len(compartmentNames_list))))
)
    Results_comp_dict={}
    for comp in particle_compartmentCoding:
        key=comp
        values=[]
        for s in Results["species"]:
            if str(particle_compartmentCoding[comp])== s[2:-7]:
                values.append(s)
            
            else:
                pass
        Results_comp_dict[key]=Results[Results["species"].isin(values)]
        
    return Results_comp_dict

def extract_by_aggSt(Results_comp_dict,MPforms_list):
    particle_forms_coding = dict(zip(MPforms_list, ["A", "B", "C", "D"]))
    Results_comp_organiced={}
    for comp in Results_comp_dict:
        Results_aggSt_dict={}
        key1=comp
        df=Results_comp_dict[comp]
        for aggst in particle_forms_coding:
            key2=aggst
            values=[]
            for c in df["species"]:
                if particle_forms_coding[aggst] == c[1]:
                    values.append(c)
                else:
                    pass
            Results_aggSt_dict[key2]=df[df["species"].isin(values)]
        Results_comp_organiced[key1]=Results_aggSt_dict
    return Results_comp_organiced    
        
def extract_by_size(Results):
    sizes=["a","b","c","d","e"]
    Results_size_dict={}
    for x in sizes:
        key=x
        values=[]
        for s in Results["species"]:
            if x in s[0:2]:
                values.append(s)
        Results_size_dict[key]= Results[Results["species"].isin(values)]
    print(Results_size_dict)
    return Results_size_dict