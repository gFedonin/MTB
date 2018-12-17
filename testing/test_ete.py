from ete3 import Tree

# we load a tree
t = Tree('((((H,K)D,(F,I)G)B,E)A,((L,(N,Q)O)J,(P,S)M)C);', format=1)
print(t)

for node in t.traverse("levelorder"):
  # Do some analysis on node
  print(node.name)

# # If we want to iterate over a tree excluding the root node, we can
# # use the iter_descendant method
# for node in t.iter_descendants("postorder"):
#   # Do some analysis on node
#   print(node.name)