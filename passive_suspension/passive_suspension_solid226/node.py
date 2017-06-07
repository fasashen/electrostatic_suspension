nodes = []
with open('./output/nodes.txt','r') as f:
        for line in f:
            node_id = int(line)
            nodes.append(node_id)
with open('./nodes_ic.txt','w') as f:
        for node in nodes:
            f.write('IC,'+str(node)+',volt,V_ic\n')