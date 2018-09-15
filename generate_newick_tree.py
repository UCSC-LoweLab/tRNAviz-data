from collections import OrderedDict, defaultdict

taxonomy = defaultdict(lambda: defaultdict(str))
nodes = []

class DAGNode(object):
  def __init__(self, name, rank, parent):
    self.name = name
    self.rank = rank
    self.parent = parent
    self.children = []

  def addChild(self, child):
    self.children.append(child)

  def isParentOf(self, node):
    return(node.name in [child.name for child in self.children])

  def indexOfChild(self, name):
    return([child.name for child in self.children].index(name))

  def __eq__(self, other):
    return(self.name == other)

  def __ne__(self, other):
    return(self.name != other)

  def __str__(self):
    children = ', '.join([child.name for child in self.children])
    return('{} ({}) -> [{}]'.format(self.name, self.rank, children))

class Tree(object):
  def __init__(self):
    self.root = DAGNode('Eukaryota', 'domain', 'root')

  # def __contains__(self, item):
  #   return(item in self)

tree = Tree()

# parse genome table into DAG
# for line in open('fungi-taxonomy-060118.tsv'):
for line in open('fungi-taxonomy-060118.tsv'):
  # build new DAG for a single species hierarchy
  assembly, taxid, varietas, species, genus, family, order, subclass, class_, subphylum, phylum, subkingdom, kingdom, dirname, dbname = line.split('\t')
  taxonomy = [('assembly', assembly),
              ('varietas', varietas),
              ('species', species),
              ('genus', genus),
              ('family', family),
              ('order', order),
              ('subclass', subclass),
              ('class', class_),
              ('subphylum', subphylum),
              ('phylum', phylum),
              ('subkingdom', subkingdom),
              ('kingdom', kingdom)]
  species_taxonomy = []
  parent = DAGNode('Eukaryota', 'domain', 'root')
  for rank, name in taxonomy[::-1]:
    if name == '': continue
    current_node = DAGNode(name, rank, parent)
    parent.addChild(current_node)
    parent = current_node


  # start from top and connect DAG to full species tree at the highest node not currently in the tree
  # initialize two pointers: current_node traverses the species DAG, while tree_node traverses the tree
  # at each iteration, check if current_node is a child of tree_node
  tree_node = tree.root
  while current_node.rank != 'domain': current_node = current_node.parent

  while True:
    if tree_node.isParentOf(current_node):
      # exit loop if terminal
      if current_node.rank == 'assembly':
        break
      # otherwise, increment both pointers
      tree_node = tree_node.children[tree_node.indexOfChild(current_node.name)]
      current_node = current_node.children[0]
      del(current_node.parent)
    else:
      tree_node.addChild(current_node)
      current_node.parent = tree_node
      break

# recursive DFS
def convertNewick(newick, node, spacer = 0):
  if len(node.children) > 0:
    children_newick = []
    for child in node.children:
      children_newick.append(convertNewick(newick, child))
    if len(children_newick) > 0:
      return('(' + ','.join(children_newick) + ')' + node.name)
  return(node.name)

print(convertNewick('', tree.root))