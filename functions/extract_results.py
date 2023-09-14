
def extract_by_size(Results):
    sizes=["a","b","c","d"]
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

def extract_by_comp(Results):
    