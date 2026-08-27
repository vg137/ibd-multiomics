import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np

#### GWAS and Parquet file iterators and converters ####

def get_gwas_iterator(gwas_filepath, chunksize=1_000_000, columns=None):
    """Returns a new iterator for the GWAS file.

    Args:
        gwas_filepath (str): The file path to the GWAS file (gzipped TSV).
        chunksize (int, optional): The number of rows to read at a time. Defaults to 10**6.
        columns (list, optional): A list of column names to read. If None, all columns are read. Defaults to None.    

    Returns:
        pandas.io.parsers.TextFileReader: An iterator that yields Pandas DataFrame chunks from the GWAS file.
    """

    gwas_file_iterator = pd.read_csv(gwas_filepath,
                                     sep='\t',
                                     usecols=columns,
                                     chunksize=chunksize,
                                    )

    return gwas_file_iterator

def get_parquet_iterator(parquet_filepath, batch_size=10**6, columns=None):
    """Returns a new iterator for the parquet file.

    Args:
        parquet_filepath (str): The file path to the Parquet file.
        batch_size (int, optional): The number of rows to read at a time. Defaults to 10**6.
        columns (list, optional): A list of column names to read. If None, all columns are read. Defaults to None.

    Returns:
        pyarrow.parquet.ParquetFile.iter_batches: An iterator that yields PyArrow RecordBatch chunks from the Parquet file.
    """

    pf = pq.ParquetFile(parquet_filepath)
    parquet_file_iterator = pf.iter_batches(columns=columns,
                                             batch_size=batch_size,
                                           )

    return parquet_file_iterator

def create_parquet_file(input_gwas_filepath, output_pf_filepath):
    """Reads data from the GWAS file in chunks and writes it to a Parquet file.

    This function iterates through chunks of the gzipped TSV GWAS file, converts each chunk
    to a PyArrow Table, and then writes these tables to a Parquet file. If the Parquet file
    does not exist, it will be created with the schema of the first chunk. Subsequent chunks
    are appended to the file.

    The output Parquet file will be compressed using Zstandard (zstd).

    Args:
        input_gwas_filepath (str): The file path to the input GWAS file (gzipped TSV).
        output_pf_filepath (str): The file path to the output Parquet file.

    Returns:
        None
    """
    writer = None

    for chunk in get_gwas_iterator(input_gwas_filepath):
        table = pa.Table.from_pandas(chunk, preserve_index=False)

        if writer is None:
            writer = pq.ParquetWriter(output_pf_filepath,
                                        table.schema,
                                        compression="zstd",
                                     )

        writer.write_table(table)

    if writer is not None:
        writer.close()

def get_snps_mask(df):
    """Returns a boolean mask indicating which rows in a DataFrame correspond to Single Nucleotide Polymorphisms (SNPs).

    A row is considered an SNP if both the 'effect_allele' and 'other_allele' values
    are single-character bases (A, C, G, or T).

    Args:
        df (pandas.DataFrame): The input DataFrame containing genetic
                               variant data, expected to have 'effect_allele' and 'other_allele' columns.

    Returns:
        pandas.Series: A boolean mask indicating which rows in the DataFrame are SNPs.
    """

    valid_bases = {"A", "C", "G", "T"}
    is_snp = (df["effect_allele"].isin(valid_bases) & df["other_allele"].isin(valid_bases))

    return is_snp

def get_num_non_snps(parquet_filepath):
    """Calculates the number of non-SNP (Single Nucleotide Polymorphism) rows in a Parquet file.

    This function iterates through batches of the Parquet file, converts each batch
    to a Pandas DataFrame, and uses `get_snps_mask` to identify and count non-SNP rows.

    Args:
        parquet_filepath (str): The file path to the Parquet file to be processed.

    Returns:
        int: The total count of non-SNP rows in the specified Parquet file.
    """

    num_non_snps = 0

    for batch in get_parquet_iterator(parquet_filepath, columns=["other_allele", "effect_allele"]):
        df = batch.to_pandas()
        df_non_snp = df[~get_snps_mask(df)]
        num_non_snps += df_non_snp.shape[0]

    return num_non_snps

def create_snp_parquet_file(input_pf_filepath, output_pf_filepath):
    """Creates a new Parquet file containing only Single Nucleotide Polymorphisms (SNPs).

    This function reads the existing Parquet file in chunks, identifies SNPs using
    the `get_snps_mask` function, and writes the filtered SNP data to a new
    Parquet file named by appending `_only_snps.parquet` to the original filename.
    The new file will have the same schema as the original, but only contain SNP rows.

    Args:
        input_pf_filepath (str): The file path to the input Parquet file.
        output_pf_filepath (str): The file path to the output Parquet file.

    Returns:
        None: The function writes the data to a file and does not return any value.
    """

    pf = pq.ParquetFile(input_pf_filepath)

    writer = None

    for batch in pf.iter_batches(batch_size=1_000_000):
        df = batch.to_pandas()

        df_snps = df[get_snps_mask(df)]

        table = pa.Table.from_pandas(df_snps, preserve_index=False)

        if writer is None:
            writer = pq.ParquetWriter(output_pf_filepath, table.schema)

        writer.write_table(table)

    if writer is not None:
        writer.close()

def get_chr_max_locations(parquet_filepath,
                            chr_colname="chromosome",
                            bp_loc_colname="base_pair_location",
                            chr_arr=range(1,24),
                         ):
    """Calculates the maximum base pair location for each chromosome in a Parquet file.

    This function iterates through a Parquet file, groups data by chromosome,
    and determines the highest base pair location for each chromosome.

    Args:
        parquet_file_path (str): The path to the Parquet file.
        chr_colname (str, optional): The name of the chromosome column. Defaults to "chromosome".
        bp_loc_colname (str, optional): The name of the base pair location column. Defaults to "base_pair_location".
        chr_arr (iterable, optional): An iterable of chromosome numbers to consider. Defaults to range(1, 24).

    Returns:
        np.ndarray: A numpy array where the index corresponds to chromosome numbers and values are their maximum base pair locations.
    """

    chr_max_loc_dict = {chr : 0 for chr in chr_arr} # dictionary to store the maximum base_pair_locations for every chromosome

    for i, batch in enumerate(get_parquet_iterator(parquet_filepath, columns=[chr_colname, bp_loc_colname])):
        for group in batch.to_pandas().groupby(chr_colname):
            chrom = group[0]
            batch_max_loc = group[1][bp_loc_colname].max()
            chr_max_loc_dict[chrom] = max(chr_max_loc_dict[chrom], batch_max_loc)

    chr_max_locs = np.array([chr_max_loc_dict[c] for c in chr_arr])

    return chr_max_locs

def get_chr_min_locations(parquet_filepath,
                            chr_colname="chromosome",
                            bp_loc_colname="base_pair_location",
                            chr_arr=range(1,24),
                         ):
    """Calculates the minimum base pair location for each chromosome in a Parquet file.

    This function iterates through a Parquet file, groups data by chromosome,
    and determines the lowest base pair location for each chromosome.

    Args:
        parquet_file_path (str): The path to the Parquet file.
        chr_colname (str, optional): The name of the chromosome column. Defaults to "chromosome".
        bp_loc_colname (str, optional): The name of the base pair location column. Defaults to "base_pair_location".
        chr_arr (iterable, optional): An iterable of chromosome numbers to consider. Defaults to range(1, 24).

    Returns:
        np.ndarray: A numpy array where the index corresponds to chromosome numbers and values are their minimum base pair locations.
    """

    chr_min_locations = {chr : np.inf for chr in chr_arr} # dictionary to store the minimum base_pair_locations for every chromosome

    for i, batch in enumerate(get_parquet_iterator(parquet_filepath, columns=[chr_colname, bp_loc_colname])):
        for group in batch.to_pandas().groupby(chr_colname):
            chrom = group[0]
            batch_min_loc = group[1][bp_loc_colname].min()
            chr_min_locations[chrom] = min(chr_min_locations[chrom], batch_min_loc)

    chr_min_locs = np.array([chr_min_locations[c] for c in chr_arr])

    return chr_min_locs

#### Gencode and SNP-to-gene mapping functions ####

def get_annotations_iterator(annotations_filepath, chunksize=100_000, columns=None):
    """Returns a new iterator for the annotation file.

    Args:
        annotations_filepath (str): The file path to the annotation file (gzipped TSV).
        chunk (int, optional): The number of rows to read at a time. Defaults to 100_000.
        columns (list, optional): A list of column names to read. If None, all columns are read. Defaults to None.

    Returns:
        pandas.io.parsers.TextFileReader: An iterator that yields Pandas DataFrame chunks from the Gencode file.
    """

    annotations_iterator = pd.read_csv(annotations_filepath,
                                        sep='\t',
                                        chunksize=chunksize,
                                        comment='#',
                                        names=["chromosome",
                                                "source",
                                                "feature",
                                                "start",
                                                "end",
                                                "score(not used)",
                                                "strand",
                                                "genomic phase",
                                                "additional_info",
                                                ],
                                        usecols=columns,
                                       )

    return annotations_iterator

def get_gene_info(additional_info_series):
    """Extracts gene information from the 'additional_info' column of a Gencode annotation DataFrame.

    This function parses the 'additional_info' strings to extract gene IDs, gene names,
    and gene types, returning them as separate Pandas Series.

    Args:
        additional_info_series (pandas.Series): A Pandas Series containing the 'additional_info' strings.
    
    Returns:
        tuple: A tuple containing three Pandas Series:
            - gene_ids (pandas.Series): Series of extracted gene IDs.
            - gene_names (pandas.Series): Series of extracted gene names.
            - gene_types (pandas.Series): Series of extracted gene types.
    """

    gene_ids = additional_info_series.str.extract(r'gene_id "([^"]+)"')[0]
    gene_names = additional_info_series.str.extract(r'gene_name "([^"]+)"')[0]
    gene_types = additional_info_series.str.extract(r'gene_type "([^"]+)"')[0]

    return gene_ids, gene_names, gene_types

def create_annotations_parquet_file(input_annotations_filepath, output_pf_filepath):
    """Reads data from the Gencode annotation file in chunks and writes it to a Parquet file.

    This function iterates through chunks of the gzipped TSV Gencode annotation file, converts each chunk
    to a PyArrow Table, and then writes these tables to a Parquet file. If the Parquet file
    does not exist, it will be created with the schema of the first chunk. Subsequent chunks
    are appended to the file.

    The output Parquet file will be compressed using Zstandard (zstd).

    Args:
        input_annotations_filepath (str): The file path to the input Gencode annotation file (gzipped TSV).
        output_pf_filepath (str): The file path to the output Parquet file.

    Returns:
        None
    """
    writer = None

    for chunk in get_annotations_iterator(input_annotations_filepath):
        chunk["gene_id"], chunk["gene_name"], chunk["gene_type"] = get_gene_info(chunk["additional_info"])
        chunk = chunk.drop(columns=["additional_info", "score(not used)", "genomic phase", "source"])
        table = pa.Table.from_pandas(chunk, preserve_index=False)

        if writer is None:
            writer = pq.ParquetWriter(output_pf_filepath,
                                        table.schema,
                                        compression="zstd",
                                     )

        writer.write_table(table)

    if writer is not None:
        writer.close()

def add_assigned_gene_info(signif_snps_df, annotations_df):
    """Assigns a gene and adds information about the assigned gene to each SNP in the given DataFrame.
    This function calculates the distance from each SNP to the transcription start site (TSS) of genes in the annotations DataFrame,
    and assigns the closest gene to each SNP. It then adds columns to the SNP DataFrame with information about the assigned gene.

    Args:
        signif_snps_df (pandas.DataFrame): A DataFrame containing SNP information with columns 'chromosome' and 'base_pair_location'.
        annotations_df (pandas.DataFrame): A DataFrame containing gene annotations with columns 'chromosome', 'tss', 'gene_id', 'gene_name', and 'gene_type'.
    
    Returns:
        None: The function modifies the signif_snps_df DataFrame in place by adding new columns.
    """

    for chrom, snp_group_df in signif_snps_df.groupby("chromosome"):
        # Filter annotations for the current chromosome
        annotations_chrom_df = annotations_df.loc[annotations_df["chromosome"] == f"chr{chrom}"]

        gene_tss_arr = annotations_chrom_df["tss"].to_numpy()
        snp_pos_arr = snp_group_df["base_pair_location"].to_numpy()

        # Calculate distances between SNPs and gene TSSs as a 2D array
        snp_tss_distances = np.abs(snp_pos_arr[:, None] - gene_tss_arr[None, :])

        closest_distance_index_arr = np.argmin(snp_tss_distances, axis=1,)
        closest_distance_arr = np.min(snp_tss_distances, axis=1,)

        assigned_gene_info = annotations_chrom_df[["gene_id", "gene_name", "gene_type"]].iloc[closest_distance_index_arr]

        signif_snps_df.loc[snp_group_df.index, ["distance_to_assigned_gene"]] = closest_distance_arr
        signif_snps_df.loc[snp_group_df.index,
                            ["assigned_gene_id", "assigned_gene_name", "assigned_gene_type"]] = assigned_gene_info.values

