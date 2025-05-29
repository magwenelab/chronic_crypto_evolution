import pandas as pd
import duckdb
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

cwd = os.getcwd()

def list_effect_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT effect_type
        FROM effects
        """
    df = con.execute(query).fetchdf()
    result = df['effect_type'].tolist()
    result.sort()
    result.insert(0, None)
    con.close()
    return result
    
def list_impacts(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT impact
        FROM effects
        """
    df = con.execute(query).fetchdf()
    result = df['impact'].tolist()
    result.sort()
    con.close()
    return result

def list_gene_names(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_name
        FROM gff
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_name'])
    result = df['gene_name'].tolist()
    result.sort()
    con.close()
    return result

def list_gene_ids(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT gene_id
        FROM gff
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['gene_id'])
    result = df['gene_id'].tolist()
    result.sort()
    con.close()
    return result

def list_samples(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT sample
        FROM samples
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['sample'])
    result = df['sample'].tolist()
    result.sort()
    con.close()
    return result

def list_strains(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT strain
        FROM samples
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['strain'])
    result = df['strain'].tolist()
    result.sort()
    con.close()
    return result

def list_lineages(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT lineage
        FROM samples
        """
    df = con.execute(query).fetchdf()
    df = df.dropna(subset=['lineage'])
    result = df['lineage'].tolist()
    result.sort()
    result.insert(0, None)
    con.close()
    return result

def effects(db,sample=None, strain=None, gene_name=None, gene_id=None, impact=None, effect_type=None, lineage=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    elif sample and strain:
        raise ValueError("Only one of Sample IDs or Strains should be provided.")
    elif (sample and lineage) or (strain and lineage):
        raise ValueError("Only one of Sample IDs/Strains or Lineage should be provided.")
    elif not (gene_name or gene_id or sample or strain or lineage):
        raise ValueError("At least one of Gene names, Gene IDs, Sample IDs, Strains or Lineage should be provided.")
    else:
        pass
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if strain:
        strain = tuple(strain)
        query_strain = f"""
            SELECT sample
            FROM samples
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_strain).fetchdf()
        sample = sample_df['sample'].tolist()
        sample = tuple(sample)
    elif sample:
        sample = tuple(sample)
    else:
        pass
    
    query = f"""
        SELECT samples.strain, presence.sample, presence.lineage,
            variants.var_id, chromosome_names.chromosome,
            variants.pos AS position, variants.ref AS reference, variants.alt AS alternative,
            effects.gene_name, effects.gene_id, effects.transcript_id,
            effects.impact, effects.effect_type, effects.effect,
            effects.codon_change, effects.amino_acid_change, effects.amino_acid_length,
            effects.transcript_biotype, effects.gene_coding, effects.exon_rank,
            mapq_depth.mean_depth_normalized, mapq_depth.mean_mapq
        FROM variants 
        JOIN chromosome_names ON variants.accession = chromosome_names.accession
        JOIN presence ON variants.var_id = presence.var_id
        JOIN effects ON variants.var_id = effects.var_id
        JOIN samples ON presence.sample = samples.sample
        LEFT JOIN mapq_depth ON mapq_depth.feature_id = effects.transcript_id AND mapq_depth.sample = presence.sample
        """
            
    if gene_name and not sample and not lineage:
        regex_pattern = '|'.join(gene_name)
        query += f"WHERE regexp_matches(effects.gene_name, '{regex_pattern}')"
    elif gene_id and not sample and not lineage:
        regex_pattern = '|'.join(gene_id)
        query += f"WHERE regexp_matches(effects.gene_id, '{regex_pattern}')"
    elif gene_name and sample:
        regex_pattern = '|'.join(gene_name)
        query += f"WHERE regexp_matches(effects.gene_name, '{regex_pattern}') AND presence.sample IN {sample}"
    elif gene_id and sample:
        regex_pattern = '|'.join(gene_id)
        query += f"WHERE regexp_matches(effects.gene_id, '{regex_pattern}') AND presence.sample IN {sample}"
    elif sample:
        query += f"WHERE presence.sample IN {sample}"
    elif (not gene_name and not gene_id) and lineage:
        lineage = tuple(lineage)
        query += f"WHERE presence.lineage IN {lineage}"
    elif gene_name and lineage:
        regex_pattern = '|'.join(gene_name)
        lineage = tuple(lineage)
        query += f"WHERE regexp_matches(effects.gene_name, '{regex_pattern}') AND presence.lineage IN {lineage}"
    elif gene_id and lineage:
        lineage = tuple(lineage)
        regex_pattern = '|'.join(gene_id)
        query += f"WHERE regexp_matches(effects.gene_id, '{regex_pattern}') AND presence.lineage IN {lineage}"
    
    if impact and effect_type:
        raise ValueError("Only one of Impacts or Effect types should be provided.")
    if impact:
        impact = tuple(impact)
        query += f"AND effects.impact IN {impact}"
    elif effect_type:
        effect_type = tuple(effect_type)
        query += f"AND effects.effect_type IN {effect_type}"
    
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

def df_to_seqrecord(df):
    records = []
    for index, row in df.iterrows():
        seq = Seq(row['seq'])
        record = SeqRecord(seq, id=f"{row['sample']}|{row['transcript_id']}", description=row['seq_description'])
        records.append(record)
    return records

def sequences(db, seq_type='DNA', sample=None, strain=None, gene_id=None, gene_name=None, filename=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    else:
        pass
    if sample and strain:
        raise ValueError("Only one of Sample IDs or Strains should be provided.")
    else:
        pass
    if not (gene_name or gene_id or sample or strain):
        raise ValueError("At least one of Gene names, Gene IDs, Sample IDs or Strains should be provided.")
    else:
        pass
    if not filename:
        raise ValueError("Filename is required.")
    else:
        pass
    
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if gene_name:
        gene_name = tuple(gene_name)
        query_gene_id = f"""
            SELECT *
            FROM gff
            WHERE gene_name IN {gene_name}
            """
        gene_id_df = con.execute(query_gene_id).fetchdf()
        gene_id = tuple(gene_id_df['gene_id'].unique().tolist())
        print(gene_id)
    elif gene_id:
        gene_id = tuple(gene_id)
        print(gene_id)

    if strain:
        strain = tuple(strain)
        query_sample = f"""
            SELECT sample
            FROM samples
            WHERE strain IN {strain}
            """
        sample_df = con.execute(query_sample).fetchdf()
        sample = tuple(sample_df['sample'].tolist())
        print(sample)
    elif sample:
        sample = tuple(sample)
        print(sample)
        
    if gene_id and not sample:
        query = f"""
            SELECT *
            FROM sequences
            WHERE transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND seq_type = '{seq_type}'"""
    elif gene_id and sample:
        query = f"""
            SELECT *
            FROM sequences
            WHERE transcript_id IN (
                SELECT DISTINCT feature_id
                FROM gff
                WHERE gene_id IN {gene_id}
            )
            AND sample IN {sample}
            AND seq_type = '{seq_type}'"""
    elif sample:
        query = f"""
            SELECT *
            FROM sequences
            WHERE sample IN {sample}
            AND seq_type = '{seq_type}'"""
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    print(f'Saving result to fasta file {filename}')
    SeqIO.write(df_to_seqrecord(result), filename, 'fasta')

def list_chromosomes(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT chromosome
        FROM chromosome_names
        """
    df = con.execute(query).fetchdf()
    result = df['chromosome'].tolist()
    result.sort()
    result.insert(0, None)
    con.close()
    return result

def get_cnv_max_length(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
    SELECT MAX("window_size") AS max_length
    FROM structural_variants    """
    df = con.execute(query).fetchdf()
    result = df['max_length'].values[0] 
    con.close()
    return result

def get_cnv(db,lineage=None, sample=None, strain=None, chromosome=None, structure=None, repeat_fraction=None, start=None, end =None,min_size=None, max_size=None):
    con = duckdb.connect(database=db, read_only=True)
    if (sample and strain) or (sample and lineage) or (strain and lineage):
        raise ValueError("Only one of Sample IDs, Strains or Lineage should be provided.")
    
    if repeat_fraction:
        q_repeat_fraction = repeat_fraction
    else:
        q_repeat_fraction = 1
    query = f"""
        SELECT chromosome_names.chromosome, structural_variants.start, structural_variants."end",
            structural_variants.structure, structural_variants.window_size, structural_variants.repeat_fraction,
            samples.sample, samples.strain, samples.lineage,
        FROM structural_variants
        JOIN samples ON structural_variants.sample = samples.sample
        JOIN chromosome_names ON structural_variants.accession = chromosome_names.accession
        WHERE structural_variants.repeat_fraction <= {q_repeat_fraction}
        """
    if lineage:
        query += f"AND samples.lineage IN {lineage}"
    if sample:
        query += f"AND samples.sample IN {sample}"
    if strain:
        query += f"AND samples.strain IN {strain}"
    if chromosome:
        query += f"AND chromosome_names.chromosome IN {chromosome} "
    if structure:
        query += f"AND structural_variants.structure == '{structure}' "
    if start:
        query += f"AND structural_variants.start >= {start} "
    if end:
        query += f"""AND structural_variants."end" <= {end} """
    if min_size:
        query += f"AND structural_variants.window_size >= {min_size} "
    if max_size:
        query += f"AND structural_variants.window_size <= {max_size} "


    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

def get_metadata(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT *
        FROM samples
        """
    df = con.execute(query).fetchdf()
    con.close()
    return df

def list_feature_types(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT primary_tag
        FROM gff
        """
    df = con.execute(query).fetchdf()
    result = df['primary_tag'].tolist()
    con.close()
    return result

def list_descriptions(db):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT DISTINCT description
        FROM gff
        """
    df = con.execute(query).fetchdf()
    result = df['description'].tolist()
    result.insert(0, None)
    con.close()
    return result

def genes(db, gene_name=None, gene_id=None, chromosome=None, start=None, end=None, feature_type=None, description=None):
    con = duckdb.connect(database=db, read_only=True)
    query = f"""
        SELECT gff.lineage, chromosome_names.chromosome,
            gff.start, gff."end", gff.primary_tag,
            gff.gene_name, gff.gene_id,
            gff.feature_id, gff.parent,
            gff.description
        FROM gff
        JOIN chromosome_names ON gff.accession = chromosome_names.accession
        """
    if gene_name:
        query += f"WHERE gene_name IN {gene_name}"
    if gene_id:
        query += f"WHERE gene_id IN {gene_id}"
    if chromosome:
        query += f"WHERE chromosome IN {chromosome}"
    if start:
        query += f"AND start >= {start}"
    if end:
        query += f"""AND "end" <= {end}"""
    if feature_type:
        query += f"AND primary_tag IN {feature_type}"
    if description:
        query += f"AND description IN {description}"
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result

def mapq_depth(db, gene_name=None, gene_id=None, feature_type=None):
    if gene_name and gene_id:
        raise ValueError("Only one of Gene names or Gene IDs should be provided.")
    elif not (gene_name or gene_id):
        raise ValueError("At least one of Gene names or Gene IDs should be provided.")
    else:
        pass
    con = duckdb.connect(database=db, read_only=True)
    con = con.execute(f"SET temp_directory = '{cwd}'")
    
    if gene_name:
        query_gene_name = f"""
            SELECT feature_id
            FROM gff
            WHERE primary_tag = 'gene' AND gene_name IN {gene_name}
            """
        feature_id_df = con.execute(query_gene_name).fetchdf()
    elif gene_id:
        query_gene_id = f"""
            SELECT feature_id
            FROM gff
            WHERE gene_id IN {gene_id} and primary_tag = 'gene'
            """
        feature_id_df = con.execute(query_gene_id).fetchdf()
    
    feature_id = tuple(feature_id_df['feature_id'].tolist())
    
    query = f"""
        SELECT *
        FROM mapq_depth
        WHERE feature_id IN {feature_id}
        """
        
    if feature_type:
        query += f"AND primary_tag IN {feature_type}"
    print(query)
    result = con.execute(query).fetchdf()
    con.close()
    return result