# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 22:53:56 2023

@author: User
"""

"""
Integrate resources to create a drug repurposing hetnet
"""

import pandas
import seaborn

import hetnetpy.hetnet #hetio package has renamed to hetnetpy
import hetnetpy.readwrite
import hetnetpy.stats

#from utils import rawgit, obo_iri #rawgit is already deprecated

import urllib.parse

def obo_iri(obo_id):
    """Construct an IRI for an OBO term ID."""
    base_url = 'http://purl.obolibrary.org/obo/'
    return urllib.parse.urljoin(base_url, obo_id.replace(':', '_'))

"""
Define the metagraph and instantiate the graph
"""

kind_to_abbev = {
    
    # metanodes
    'Compound': 'C',
    'Disease': 'D',
    'Gene': 'G',
    'Anatomy': 'A',
    'Symptom': 'S',
    'Side Effect': 'SE',
    'Pathway': 'PW',
    'Pharmacologic Class': 'PC',
    'Biological Process': 'BP',
    'Cellular Component': 'CC',
    'Molecular Function': 'MF',
    
    # metaedges
    'treats': 't',
    'palliates': 'p',
    'binds': 'b',
    'expresses': 'e',
    'regulates': 'r',
    'upregulates': 'u',
    'downregulates': 'd',
    'interacts': 'i',
    'includes': 'i',
    'covaries': 'c',
    'regulates': 'r',
    'participates': 'p',
    'resembles': 'r',
    'associates': 'a',
    'localizes': 'l',
    'presents': 'p',
    'causes': 'c',
}

metaedge_tuples = [
    ('Compound', 'Disease', 'treats', 'both'),
    ('Compound', 'Disease', 'palliates', 'both'),
    ('Compound', 'Gene', 'binds', 'both'),
    ('Compound', 'Gene', 'upregulates', 'both'),
    ('Compound', 'Gene', 'downregulates', 'both'),
    ('Compound', 'Compound', 'resembles', 'both'),
    ('Compound', 'Side Effect', 'causes', 'both'),
    ('Pharmacologic Class', 'Compound', 'includes', 'both'),
    ('Anatomy', 'Gene', 'expresses', 'both'),
    ('Anatomy', 'Gene', 'upregulates', 'both'),
    ('Anatomy', 'Gene', 'downregulates', 'both'),
    ('Gene', 'Gene', 'interacts', 'both'),
    ('Gene', 'Gene', 'covaries', 'both'),
    ('Gene', 'Gene', 'regulates', 'forward'),
    ('Gene', 'Pathway', 'participates', 'both'),
    ('Gene', 'Biological Process', 'participates', 'both'),
    ('Gene', 'Cellular Component', 'participates', 'both'),
    ('Gene', 'Molecular Function', 'participates', 'both'),
    ('Disease', 'Disease', 'resembles', 'both'),
    ('Disease', 'Gene', 'associates', 'both'),
    ('Disease', 'Gene', 'upregulates', 'both'),
    ('Disease', 'Gene', 'downregulates', 'both'),
    ('Disease', 'Anatomy', 'localizes', 'both'),
    ('Disease', 'Symptom', 'presents', 'both'),
]
metagraph = hetnetpy.hetnet.MetaGraph.from_edge_tuples(metaedge_tuples, kind_to_abbev)
graph = hetnetpy.hetnet.Graph(metagraph)

""" 
Gene Nodes
"""

commit = 'a7362748a34211e5df6f2d185bb3246279760546'
owner = 'dhimmel'
repo = 'entrez-gene'
path = 'data'
filename = 'genes-human.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

gene_df = pandas.read_table(url)
gene_df = gene_df[gene_df.type_of_gene == 'protein-coding']
coding_genes = set(gene_df.GeneID)
gene_df.head(2)

for i, row in gene_df.iterrows():
    if row.type_of_gene != 'protein-coding':
        continue
    data = {
        'description': row['description'],
        'source': 'Entrez Gene',
        'url': 'http://identifiers.org/ncbigene/{}'.format(row.GeneID),
        'license': 'CC0 1.0',
    }
    if pandas.notnull(row['chromosome']):
        data['chromosome'] = row['chromosome']
    graph.add_node(kind = 'Gene', identifier=row.GeneID, name=row.Symbol, data=data)
    
"""
Diseases Nodes
"""

commit = '75050ea2d4f60e745d3f3578ae03560a2cc0e444'
owner = 'dhimmel'
repo = 'disease-ontology'
path = 'data'
filename = 'slim-terms.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

disease_df = pandas.read_table(url)
disease_df.head(2)

for i, row in disease_df.iterrows():
    data = {
        'source': 'Disease Ontology',
        'url': obo_iri(row.doid),
        'license': 'CC BY 3.0',
    }
    graph.add_node(kind='Disease', identifier=row.doid, name=row['name'], data=data)
    

"""
Compound Nodes
"""
commit = '3e87872db5fca5ac427ce27464ab945c0ceb4ec6'
owner = 'dhimmel'
repo = 'drugbank'
path = 'data'
filename = 'drugbank-slim.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

compound_df = pandas.read_table(url)
compound_df.head(2)

for i, row in compound_df.iterrows():
    url = 'http://www.drugbank.ca/drugs/' + row.drugbank_id
    data = {
        'source': 'DrugBank',
        'inchikey': row.inchikey,
        'inchi': row.inchi, 'url': url,
        'license': 'CC BY-NC 4.0',
    }
    graph.add_node(kind='Compound', identifier=row.drugbank_id, name=row['name'], data=data)
    
"""
Anotomy Nodes
"""

commit = '134f23479186abba03ba340fc6dc90e16c781920'
owner = 'dhimmel'
repo = 'uberon'
path = 'data'
filename = 'hetio-slim.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

uberon_df = pandas.read_table(url)
uberon_df.head(2)

for i, row in uberon_df.iterrows():
    data = {
        'source': 'Uberon',
        'url': obo_iri(row['uberon_id']),
        'license': 'CC BY 3.0',
    }
    for xref in 'mesh_id', 'bto_id':
        if pandas.notnull(row[xref]):
            data[xref] = row[xref]
    graph.add_node(kind='Anatomy', identifier=row['uberon_id'], name=row['uberon_name'], data=data)
    
""" 
Symptom Nodes
"""

commit = 'a7036a37302973b15ab949aab4056d9bc062910e'
owner = 'dhimmel'
repo = 'mesh'
path = 'data'
filename = 'symptoms.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

symptom_df = pandas.read_table(url)
symptom_df.head(2)

for i, row in symptom_df.iterrows():
    url = 'http://identifiers.org/mesh/{}'.format(row.mesh_id)
    data = {
        'source': 'MeSH',
        'url': url,
        'license': 'CC0 1.0',
    }
    graph.add_node(kind='Symptom', identifier=row.mesh_id, name=row.mesh_name, data=data)
    
"""
Pathway Nodes and Edges
"""

commit = '1bd2c68853e38297d20f8f885419ff81fc0608a8'
owner = 'dhimmel'
repo = 'pathways'
path = 'data'
filename = 'pathways.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

pathway_df = pandas.read_table(url)
pathway_df = pathway_df.query("n_coding_genes > 1")
source_map = {
    'wikipathways': 'WikiPathways',
    'reactome': 'Reactome via Pathway Commons',
    'pid': 'PID via Pathway Commons',
}
pathway_df.source = pathway_df.source.map(source_map)
pathway_df.tail(2)

for i, row in pathway_df.iterrows():
    pathway_id = row.identifier
    data = {'license': row.license, 'source': row.source}
    if pandas.notnull(row.url):
        data['url'] = row.url
    graph.add_node(kind='Pathway', identifier=pathway_id, name=row['name'], data=data)
    for gene in row.coding_genes.split('|'):
        gene = int(gene)
        source_id = 'Gene', gene
        target_id = 'Pathway', pathway_id
        edge_data = data.copy()
        edge_data['unbiased'] = False
        graph.add_edge(source_id, target_id, 'participates', 'both', edge_data)
        
"""
Pharmacologic Classes
"""

commit = 'e80a0c966a53ce48650d98069b126801c2793517'
owner = 'dhimmel'
repo = 'drugcentral'
path = 'rephetio'
filename = 'classes.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

class_df = pandas.read_table(url)
class_types = {'Physiologic Effect', 'Mechanism of Action', 'Chemical/Ingredient'}
class_df = class_df.query("class_type in @class_types")
class_df.head(2)

for row in class_df.itertuples():
    data = {
        'source': '{} via DrugCentral'.format(row.class_source),
        'class_type': row.class_type,
        'license': 'CC BY 4.0',
        'url': row.url,
    }
    graph.add_node(kind='Pharmacologic Class', identifier=row.class_id, name=row.class_name, data=data)

filename = 'drug-to-class.tsv'
url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'
drug_class_df = pandas.read_table(url)
drug_class_df = drug_class_df.query("class_id in @class_df.class_id")
drug_class_df.head(2)

for row in drug_class_df.itertuples():
    data = {
        'source': 'DrugCentral',
        'license': 'CC BY 4.0',
        'unbiased': False,
    }
    source_id = 'Pharmacologic Class', row.class_id
    target_id = 'Compound', row.drugbank_id
    graph.add_edge(source_id, target_id, 'includes', 'both', data)
### Assertion Error: edge already exist 

""" 
Gen Ontology Domains
"""

commit = '87bab297f55db283e65a7a984607316b409415ae'
owner = 'dhimmel'
repo = 'gene-ontology'
path = 'annotations/taxid_9606'
filename = 'GO_annotations-9606-inferred-allev.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

go_df = pandas.read_table(url)
go_df.head(2)

for i, row in go_df.iterrows():
    genes = coding_genes & set(map(int, row.gene_ids.split('|')))
    if 2 > len(genes) or len(genes) > 1000:
        continue
    kind = row['go_domain'].replace('_', ' ').title()
    data = {'source': 'Gene Ontology', 'url': obo_iri(row.go_id), 'license': 'CC BY 4.0'}
    target = graph.add_node(kind=kind, identifier=row['go_id'], name=row['go_name'], data=data)
    target_id = target.get_id()
    for gene in genes:
        source_id = 'Gene', gene
        data = {'source': 'NCBI gene2go', 'unbiased': False, 'license': 'CC BY 4.0'}
        graph.add_edge(source_id, target_id, 'participates', 'both', data)
        
"""
Disease-gene associations from compilation
"""

commit = '3141b25edf91747accfece5e297c77ab1b37f4da'
owner = 'dhimmel'
repo = 'integrate'
path = 'compile'
filename = 'DaG-association.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

association_df = pandas.read_table(url)
association_df.head(2)

for i, row in association_df.iterrows():
    source_id = 'Disease', row.doid_id
    target_id = 'Gene', row.entrez_gene_id
    sources = sorted(row.sources.split('|'))
    data = {'sources': sources, 'unbiased': 'GWAS Catalog' in sources}
    if pandas.notnull(row['license']):
        data['license'] = row['license']
    graph.add_edge(source_id, target_id, 'associates', 'both', data)
### KeyError: ('Disease', 'DOID:2531')
    
"""
Disease-gene differential expression
"""

commit = '1a11633b5e0095454453335be82012a9f0f482e4'
owner = 'dhimmel'
repo = 'stargeo'
path = 'data'
filename = 'diffex.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

stargeo_df = pandas.read_table(url)
# Filter to at most 250 up and 250 down-regulated genes per disease
stargeo_df = stargeo_df.groupby(['slim_id', 'direction']).apply(
    lambda df: df.nsmallest(250, 'p_adjusted')).reset_index(drop=True)
stargeo_df.head(2)

for row in stargeo_df.itertuples():
    source_id = 'Disease', row.slim_id
    target_id = 'Gene', row.entrez_gene_id
    kind = row.direction + 'regulates'
    data = {
        'source': 'STARGEO',
        'log2_fold_change': round(row.log2_fold_change, 5),
        'unbiased': True,
        'license': 'CC0 1.0'
    }
    graph.add_edge(source_id, target_id, kind, 'both', data)
### KeyError: ('Disease', 'DOID:0050156')

"""
Chemical similarity
"""

commit = '3e87872db5fca5ac427ce27464ab945c0ceb4ec6'
owner = 'dhimmel'
repo = 'drugbank'
path = 'data'
filename = 'similarity-slim.tsv.gz'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

chemical_df = pandas.read_table(url, compression='gzip')
chemical_df = chemical_df[chemical_df.similarity >= 0.5]
chemical_df.head(2)

for i, row in chemical_df.iterrows():
    source_id = 'Compound', row.compound0
    target_id = 'Compound', row.compound1
    data = {
        'source': 'Dice similarity of ECFPs',
        'similarity': round(row.similarity, 4),
        'unbiased': True,
        'license': 'CC0 1.0',
    }
    graph.add_edge(source_id, target_id, 'resembles', 'both', data)
    
"""
MEDLINE cooccurrence data
"""

medline_data = {
    'source': 'MEDLINE cooccurrence',
    'unbiased': False,
    'license': 'CC0 1.0',
}

"""
Symptom edges
"""

commit = '60d611892bf387b5b23c5f2e2e3bc472cfce85f3'
owner = 'dhimmel'
repo = 'medline'
path = 'data'
filename = 'disease-symptom-cooccurrence.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

disease_symptom_df = pandas.read_table(url)
disease_symptom_df = disease_symptom_df[disease_symptom_df.p_fisher < 0.005]
disease_symptom_df.head(2)

for i, row in disease_symptom_df.iterrows():
    source_id = 'Disease', row.doid_code
    target_id = 'Symptom', row.mesh_id
    data = medline_data.copy()
    graph.add_edge(source_id, target_id, 'presents', 'both', data)
### KeyError: ('Disease', 'DOID:10652')

"""
Disease-localization edges
"""

commit = '60d611892bf387b5b23c5f2e2e3bc472cfce85f3'
owner = 'dhimmel'
repo = 'medline'
path = 'data'
filename = 'disease-uberon-cooccurrence.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

disease_anatomy_df = pandas.read_table(url)
disease_anatomy_df = disease_anatomy_df[disease_anatomy_df.p_fisher < 0.005]
disease_anatomy_df = disease_anatomy_df[disease_anatomy_df.uberon_id.isin(uberon_df['uberon_id'])]
disease_anatomy_df.head(2)

for i, row in disease_anatomy_df.iterrows():
    source_id = 'Disease', row.doid_code
    target_id = 'Anatomy', row.uberon_id
    data = medline_data.copy()
    graph.add_edge(source_id, target_id, 'localizes', 'both', data)
### KeyError: ('Disease', 'DOID:10652')

"""
Disease-disease similarity
"""

commit = '60d611892bf387b5b23c5f2e2e3bc472cfce85f3'
owner = 'dhimmel'
repo = 'medline'
path = 'data'
filename = 'disease-disease-cooccurrence.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

disease_similarity_df = pandas.read_table(url)
disease_similarity_df = disease_similarity_df[-disease_similarity_df[['doid_code_0', 'doid_code_1']].apply(frozenset, 1).duplicated()]
disease_similarity_df = disease_similarity_df[disease_similarity_df.p_fisher < 0.005]
disease_similarity_df.head(2)

for i, row in disease_similarity_df.iterrows():
    source_id = 'Disease', row.doid_code_0
    target_id = 'Disease', row.doid_code_1
    data = medline_data.copy()
    graph.add_edge(source_id, target_id, 'resembles', 'both', data)
### KeyError: ('Disease', 'DOID:10652')

"""
Anatomy-gene expression presence
"""

commit = '3141b25edf91747accfece5e297c77ab1b37f4da'
owner = 'dhimmel'
repo = 'integrate'
path = 'compile'
filename = 'AeG-expression.tsv.gz'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

expr_df = pandas.read_table(url, low_memory=False)
expr_df = expr_df[expr_df.uberon_id.isin(uberon_df.uberon_id) & expr_df.entrez_gene_id.isin(coding_genes)]
expr_df.head(2)

for i, row in expr_df.iterrows():
    source_id = 'Gene', row['entrez_gene_id']
    target_id = 'Anatomy', row['uberon_id']
    data = {'unbiased': bool(row['unbiased'])}
    if pandas.notnull(row['license']):
        data['license'] = row['license']
    data['sources'] = row['sources'].split('|')
    graph.add_edge(source_id, target_id, 'expresses', 'both', data)

"""
Anatomy-gene differential expression
"""

commit = '08ba54e83ee8e28dec22b4351d29e23f1d034d30'
owner = 'dhimmel'
repo = 'bgee'
path = 'data'
filename = 'diffex.tsv.gz'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

diffex_df = pandas.read_table(url, compression='gzip')
diffex_df = pandas.melt(diffex_df, id_vars='GeneID', var_name='uberon_id', value_name='direction')
diffex_df = diffex_df.query("direction != 0")
diffex_df = diffex_df[diffex_df.uberon_id.isin(uberon_df.uberon_id) & diffex_df.GeneID.isin(coding_genes)]
diffex_df = diffex_df.replace({'direction': {-1: 'downregulates', 1: 'upregulates'}})
diffex_df.head(2)

for i, row in diffex_df.iterrows():
    source_id = 'Gene', row['GeneID']
    target_id = 'Anatomy', row['uberon_id']
    data = {'source': 'Bgee', 'unbiased': True}
    graph.add_edge(source_id, target_id, row['direction'], 'both', data)
    
"""
Compound bindings
"""

commit = '3141b25edf91747accfece5e297c77ab1b37f4da'
owner = 'dhimmel'
repo = 'integrate'
path = 'compile'
filename = 'CbG-binding.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

binding_df = pandas.read_table(url)
binding_df = binding_df.merge(compound_df[['drugbank_id']])
binding_df = binding_df[binding_df.entrez_gene_id.isin(coding_genes)]
binding_df.head(2)

for i, row in binding_df.iterrows():
    source_id = 'Compound', row.drugbank_id
    target_id = 'Gene', row.entrez_gene_id
    data = {'unbiased': False}
    # singular fields
    for key in 'affinity_nM', 'license':
        value = row[key]
        if pandas.notnull(value):
            data[key] = value
    # compound fields
    for key in 'sources', 'pubmed_ids', 'actions', 'urls':
        value = row[key]
        if pandas.notnull(value):
            data[key] = value.split('|')
    graph.add_edge(source_id, target_id, 'binds', 'both', data)
    
"""
Protein Interactions
"""

commit = 'f6a7edbc8de6ba2d7fe1ef3fee4d89e5b8d0b900'
owner = 'dhimmel'
repo = 'ppi'
path = 'data'
filename = 'ppi-hetio-ind.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

ppi_df = pandas.read_table(url)
ppi_df = ppi_df[ppi_df.gene_0.isin(coding_genes) & ppi_df.gene_1.isin(coding_genes)]
ppi_df.head(2)

for i, row in ppi_df.iterrows():
    source_id = 'Gene', row.gene_0
    target_id = 'Gene', row.gene_1
    data = {
        'sources': row.sources.split('|'),
        'unbiased': bool(row.unbiased),
    }
    graph.add_edge(source_id, target_id, 'interacts', 'both', data)

"""
Evolutionary rate covariation
"""

commit = '757733f77a89499439c887acb88456e011c5322e'
owner = 'dhimmel'
repo = 'erc'
path = 'data'
filename = 'erc_mam33-entrez-gt-0.6.tsv.gz'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

erc_df = pandas.read_table(url, compression='gzip')
erc_df = erc_df[erc_df.correlation >= 0.75]
erc_df = erc_df[erc_df.source_entrez.isin(coding_genes) & erc_df.target_entrez.isin(coding_genes)]
erc_df.head(2)

for i, row in erc_df.iterrows():
    source_id = 'Gene', row.source_entrez
    target_id = 'Gene', row.target_entrez
    data = {
        'source': 'ERC',
        'unbiased': True,
    }
    graph.add_edge(source_id, target_id, 'covaries', 'both', data)
### AssertionError: edge already exists

"""
Indications from the PharmacotherapyDB
"""
commit = '11d535ba0884ee56c3cd5756fdfb4985f313bd80'
owner = 'dhimmel'
repo = 'indications'
path = 'catalog'
filename = 'indications.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

indication_df = pandas.read_table(url)
categories = {'DM', 'SYM'}
indication_df = indication_df.query("category in @categories").copy()
indication_df['kind'] = indication_df.category.map({'DM': 'treats', 'SYM': 'palliates'})
indication_df.head(2)

for i, row in indication_df.iterrows():
    source_id = 'Disease', row.doid_id
    target_id = 'Compound', row.drugbank_id
    data = {'source': 'PharmacotherapyDB', 'unbiased': False, 'license': 'CC0 1.0'}
    graph.add_edge(source_id, target_id, row['kind'], 'both', data)
### KeyError: ('Disease', 'DOID:10652')

"""
LINCS L1000 relationships
"""

commit = 'abcb12f942f93e3ee839e5e3593f930df2c56845'
def filter_l1000_df(df, n):
    """
    Filter LINCS L1000 differentially expression genes to at most `n` genes
    per perturbagen-direction-status combination.
    """
    df = df.groupby(['perturbagen', 'direction', 'status']).apply(
        lambda x: x.nlargest(n, 'nlog10_bonferroni_pval')).reset_index(drop=True)
    return df

"""
LINCS compound-gene dysregulation
"""

owner = 'dhimmel'
repo = 'lincs'
path = 'data/consensi/signif'
filename = 'dysreg-drugbank.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

l1000_df = pandas.read_table(url)
l1000_df = l1000_df.query("perturbagen in @compound_df.drugbank_id and entrez_gene_id in @coding_genes")
l1000_df = filter_l1000_df(l1000_df, n=125)
l1000_df.tail(2)

mapper = {'up': 'upregulates', 'down': 'downregulates'}
for row in l1000_df.itertuples():
    source_id = 'Compound', row.perturbagen
    target_id = 'Gene', row.entrez_gene_id
    data = {
        'source': 'LINCS L1000',
        'z_score': round(row.z_score, 3),
        'method': row.status,
        'unbiased': True,
    }
    kind = mapper[row.direction]
    graph.add_edge(source_id, target_id, kind, 'both', data)
    
"""
LINCS genetic perturbations
"""

owner = 'dhimmel'
repo = 'lincs'
path = 'data/consensi/signif/'
filename = 'dysreg-knockdown.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

l1000_kd_df = filter_l1000_df(pandas.read_table(url), n=50)

owner = 'dhimmel'
repo = 'lincs'
path = 'data/consensi/signif/'
filename = 'dysreg-overexpression.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

l1000_oe_df = filter_l1000_df(pandas.read_table(url), n=50)

mapper = {'up': 'knockdown upregulates', 'down': 'knockdown downregulates'}
l1000_kd_df['kind'] = l1000_kd_df.direction.map(lambda x: mapper[x])

mapper = {'up': 'overexpression upregulates', 'down': 'overexpression downregulates'}
l1000_oe_df['kind'] = l1000_oe_df.direction.map(lambda x: mapper[x])

l1000_genetic_df = pandas.concat([l1000_kd_df, l1000_oe_df])
l1000_genetic_df = l1000_genetic_df.query('perturbagen in @coding_genes and entrez_gene_id in @coding_genes')
l1000_genetic_df = l1000_genetic_df.query('perturbagen != entrez_gene_id')

l1000_genetic_df.head(2)

for (pert, gene), df in l1000_genetic_df.groupby(['perturbagen', 'entrez_gene_id'], sort=False):
    source_id = 'Gene', pert
    target_id = 'Gene', gene
    method, = df.status.unique()
    data = {'source': 'LINCS L1000', 'subtypes': list(df.kind), 'method': method, 'unbiased': True}
    graph.add_edge(source_id, target_id, 'regulates', 'forward', data)
    
"""
Side Effects - SIDER
"""

commit = 'be3adebc0d845baaddb907a880890cb5e85f5801'
owner = 'dhimmel'
repo = 'SIDER4'
path = 'data'
filename = 'side-effect-terms.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

side_effect_df = pandas.read_table(url)
side_effect_df.head(2)

for i, row in side_effect_df.iterrows():
    umls_id = row['umls_cui_from_meddra']
    data = {
        'source': 'UMLS via SIDER 4.1',
        'url': 'http://identifiers.org/umls/{}'.format(umls_id),
        'license': 'CC BY-NC-SA 4.0',
    }
    graph.add_node(kind='Side Effect', identifier=umls_id, name=row['side_effect_name'], data=data)

owner = 'dhimmel'
repo = 'SIDER4'
path = 'data'
filename = 'side-effects.tsv'

url = f'https://raw.githubusercontent.com/{owner}/{repo}/{commit}/{path}/{filename}'

sider_df = pandas.read_table(url)
sider_df = sider_df[sider_df.drugbank_id.isin(compound_df.drugbank_id)]
sider_df.head(2)

for i, row in sider_df.iterrows():
    umls_id = row.umls_cui_from_meddra
    source_id = 'Compound', row.drugbank_id
    target_id = 'Side Effect', umls_id
    data = {
        'source': 'SIDER 4.1',
        'url': 'http://sideeffects.embl.de/se/{}/'.format(umls_id),
        'unbiased': False,
        'license': 'CC BY-NC-SA 4.0',
    }
    graph.add_edge(source_id, target_id, 'causes', 'both', data)
    
"""
Network visualizations and stats
"""

# Export node degree tables
hetnetpy.stats.degrees_to_excel(graph, 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/degrees.xlsx')

# Create and save degree distribution vizualizations
hetnetpy.stats.plot_degrees(graph, 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/degrees.pdf')
### ValueError: Number of rows must be a positive integer, not 0

# Summary of metanodes and cooresponding nodes
metanode_df = hetnetpy.stats.get_metanode_df(graph)
metanode_df.to_csv('C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/metanodes.tsv', sep='\t', index=False)
metanode_df

# Total number of nodes
metanode_df.nodes.sum()

# Summary of metaedges and cooresponding edges
metaedge_df = hetnetpy.stats.get_metaedge_df(graph)

# Calculate number of unbiased edges
rows = list()
for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
    unbiased = sum(edge.data['unbiased'] for edge in edges)
    rows.append({'metaedge': str(metaedge), 'unbiased': unbiased})

metaedge_df = metaedge_df.merge(pandas.DataFrame(rows))
metaedge_df.to_csv('C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/metaedges.tsv', sep='\t', index=False)
metaedge_df

# Summary of different styles for representing each metaedge
metaedge_style_df = hetnetpy.stats.get_metaedge_style_df(metagraph)
metaedge_style_df.to_csv('C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/metaedge-styles.tsv', sep='\t', index=False)

# Number of edges in the network
metaedge_df.edges.sum()

# Write nodes to a table
path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/nodes.tsv'
hetnetpy.readwrite.write_nodetable(graph, path)
### UnicodeEncodeError: 'charmap' codec can't encode character '\u03b1' in position 48: character maps to <undefined>

# Write edges to a table
path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/edges.sif.gz'
hetnetpy.readwrite.write_sif(graph, path)

# Write a subset of edges to a table
path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/edges-10.sif'
hetnetpy.readwrite.write_sif(graph, path, max_edges=10)

path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/edges-5k.sif.gz'
hetnetpy.readwrite.write_sif(graph, path, max_edges=5000)

# Write metagraph as json
path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/metagraph.json'
hetnetpy.readwrite.write_metagraph(metagraph, path)

# Write graph as json
path = 'C:/Users/User/OneDrive - National University of Singapore/Capstone_KG/hetnet.json.bz2'
hetnetpy.readwrite.write_graph(graph, path)

# """! sha256sum data/hetnet.json.bz2""" don know what is this for 


""" 
Barplots of metaedge and metanode counts
"""

%matplotlib inline

ax = seaborn.barplot(x='metanode', y='nodes', data=metanode_df.sort_values('nodes'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('nodes');

ax = seaborn.barplot(x='metaedge', y='edges', data=metaedge_df.sort_values('edges'))
for tick in ax.get_xticklabels():
    tick.set_rotation(90)
ax.set_xlabel(''); ax.set_ylabel('edges');

