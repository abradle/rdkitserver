
def remove_keys(in_list): 
    """Prepares a dict to be returned as JSON"""
    banned_keys = ["RDMOL", "FP"]
    out_list = []
    for item in in_list:
        for key in banned_keys:
            if key in item:
                del item[key]
        out_list.append(item)
    return out_list
 
