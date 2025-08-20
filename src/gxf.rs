//! The fastest GTF/GFF-to-BED converter chilling around
//! Alejandro Gonzales-Irribarren, 2025

mod attr;
pub use attr::*;

use std::collections::BTreeSet;

/// Represents a record parsed from a GXF (GTF/GFF) file.
///
/// This struct holds the essential fields of a GXF record, including chromosomal location,
/// feature type, genomic coordinates, strand, frame, and a parsed attribute map.
///
/// # Fields
/// * `chr` - The name of the chromosome or sequence contig.
/// * `feature` - The feature type (e.g., "gene", "transcript", "exon").
/// * `start` - The 0-based start coordinate of the feature.
/// * `end` - The 1-based end coordinate of the feature.
/// * `strand` - The strand of the feature (Forward, Reverse, or Unknown).
/// * `frame` - The frame of the feature (e.g., "0", "1", "2", ".").
/// * `attr` - A parsed `Attribute` map containing key-value pairs from the last field.
#[derive(Debug, PartialEq)]
pub struct GxfRecord<'a> {
    pub chr: String,
    pub feature: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub frame: String,
    pub attr: Attribute<'a>,
}

/// Represents the strand of a genomic feature.
///
/// This enum indicates whether a feature is on the forward (+), reverse (-), or an unknown (.) strand.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl std::fmt::Display for Strand {
    /// Formats the `Strand` enum into its string representation.
    ///
    /// `Strand::Forward` becomes `"+"`, `Strand::Reverse` becomes `"-"`, and `Strand::Unknown` becomes `"."`.
    ///
    /// # Arguments
    /// * `f` - The formatter.
    ///
    /// # Returns
    /// A `std::fmt::Result` indicating success or failure.
    ///
    /// # Example
    /// ```rust
    /// use gxf2bed::Strand; //  with the actual crate name
    /// assert_eq!(format!("{}", Strand::Forward), "+");
    /// assert_eq!(format!("{}", Strand::Reverse), "-");
    /// assert_eq!(format!("{}", Strand::Unknown), ".");
    /// ```
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

impl std::str::FromStr for Strand {
    type Err = &'static str;

    /// Parses a string slice into a `Strand` enum.
    ///
    /// Recognizes `"+"` for `Forward`, `"-"` for `Reverse`, and any other string
    /// (including an empty string) as `Unknown`.
    ///
    /// # Arguments
    /// * `s` - The string slice to parse.
    ///
    /// # Returns
    /// A `Result` containing the parsed `Strand` or a static string error if parsing fails.
    ///
    /// # Example
    /// ```rust
    /// use std::str::FromStr;
    /// use gxf2bed::Strand; //  with the actual crate name
    ///
    /// assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
    /// assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
    /// assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);
    /// assert_eq!(Strand::from_str("anything_else").unwrap(), Strand::Unknown);
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            _ => Ok(Strand::Unknown),
        }
    }
}

/// Represents the type of a gene prediction record.
///
/// This enum helps categorize gene prediction records as either parent, child, or unknown,
/// useful in hierarchical gene structures.
///
/// # Variants
/// * `Parent` - Indicates a parent record (e.g., a gene).
/// * `Child` - Indicates a child record (e.g., an exon or transcript belonging to a gene).
/// * `Unknown` - The record type is not specified or recognized.
#[derive(Debug, PartialEq)]
pub enum RecordType {
    Parent,
    Child,
    Unknown,
}

/// Represents a gene prediction, often derived from GXF records.
///
/// This struct aggregates information about a gene or transcript, including its genomic
/// coordinates, strand, and a collection of its exons. It also includes a `record_type`
/// to distinguish between parent (e.g., gene) and child (e.g., transcript) entries.
///
/// # Fields
/// * `chr` - The chromosome name.
/// * `start` - The 0-based start coordinate of the gene prediction.
/// * `end` - The 1-based end coordinate of the gene prediction.
/// * `strand` - The strand of the gene prediction.
/// * `exons` - A `BTreeSet` of (start, end) tuples representing the exons.
/// * `record_type` - The type of record (Parent, Child, or Unknown).
/// * `blocks` - An optional `BTreeSet` of (start, end) tuples for additional blocks.
#[derive(Debug, PartialEq)]
pub struct GenePred {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: BTreeSet<(u64, u64)>,
    pub record_type: RecordType,
    pub blocks: Option<BTreeSet<(u64, u64)>>,
}

impl GenePred {
    /// Creates a new, empty `GenePred` instance with default values.
    ///
    /// # Returns
    /// A new `GenePred` instance.
    ///
    /// # Example
    /// ```rust
    /// use gxf2bed::{GenePred, Strand, RecordType}; //
    /// let gene_pred = GenePred::new();
    /// assert_eq!(gene_pred.chr, "");
    /// assert_eq!(gene_pred.start, 0);
    /// assert_eq!(gene_pred.strand, Strand::Unknown);
    /// assert!(gene_pred.exons.is_empty());
    /// assert_eq!(gene_pred.record_type, RecordType::Unknown);
    /// ```
    pub fn new() -> Self {
        Self {
            chr: String::new(),
            start: 0,
            end: 0,
            strand: Strand::Unknown,
            exons: BTreeSet::new(),
            record_type: RecordType::Unknown,
            blocks: None,
        }
    }

    /// Merges another `GenePred` instance into the current one.
    ///
    /// This method updates the `GenePred` based on the `record_type` of the `query` instance.
    ///
    /// - If `query` is `RecordType::Parent`, it overwrites chromosome, start, end, strand,
    ///   and record type, and extends the exons.
    /// - If `query` is `RecordType::Child` or `RecordType::Unknown`, it updates the
    ///   chromosome, strand, and end if they are initially empty or smaller,
    ///   adjusts the start coordinate to be the minimum, and extends the exons.
    ///   The `record_type` is set to `Child` if it wasn't already `Parent`.
    ///
    /// # Arguments
    /// * `query` - The `GenePred` instance to merge from.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType}; //
    ///
    /// let mut gp1 = GenePred::new();
    /// gp1.chr = "chr1".to_string();
    /// gp1.start = 100;
    /// gp1.end = 200;
    /// gp1.strand = Strand::Forward;
    /// gp1.exons.insert((100, 150));
    /// gp1.record_type = RecordType::Child;
    ///
    /// let mut gp2 = GenePred::new();
    /// gp2.chr = "chr1".to_string();
    /// gp2.start = 50;
    /// gp2.end = 250;
    /// gp2.strand = Strand::Forward;
    /// gp2.exons.insert((50, 80));
    /// gp2.exons.insert((200, 250));
    /// gp2.record_type = RecordType::Child;
    ///
    /// gp1.merge(gp2);
    ///
    /// assert_eq!(gp1.chr, "chr1");
    /// assert_eq!(gp1.start, 50);
    /// assert_eq!(gp1.end, 250);
    /// assert_eq!(gp1.record_type, RecordType::Child);
    /// let mut expected_exons = BTreeSet::new();
    /// expected_exons.insert((50, 80));
    /// expected_exons.insert((100, 150));
    /// expected_exons.insert((200, 250));
    /// assert_eq!(gp1.exons, expected_exons);
    /// ```
    pub fn merge(&mut self, query: GenePred) {
        match query.record_type {
            RecordType::Parent => {
                self.chr = query.chr;
                self.start = query.start;
                self.end = query.end;
                self.strand = query.strand;
                self.record_type = RecordType::Parent;
                self.exons.extend(query.exons);
            }
            RecordType::Child => {
                // INFO: covers empty cases in reducing step!
                if self.chr.is_empty() {
                    self.chr = query.chr;
                    self.strand = query.strand;
                    self.end = self.end.max(query.end);
                }

                if self.start < 1 && query.start > 0 {
                    self.start = query.start;
                } else if self.start > 0 && query.start > 0 {
                    self.start = self.start.min(query.start);
                }

                self.exons.extend(query.exons);
                if self.record_type != RecordType::Parent {
                    self.record_type = RecordType::Child;
                }
            }
            RecordType::Unknown => {
                if self.chr.is_empty() && !query.chr.is_empty() {
                    self.chr = query.chr;
                    self.strand = query.strand;
                    self.end = self.end.max(query.end);
                }

                if self.start < 1 && query.start > 0 {
                    self.start = query.start;
                } else if self.start > 0 && query.start > 0 {
                    self.start = self.start.min(query.start);
                }
            }
        }
    }

    /// Returns the number of exons in the gene prediction.
    ///
    /// # Returns
    /// The count of exons as `usize`.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType}; //
    /// let mut gp = GenePred::new();
    /// gp.exons.insert((100, 200));
    /// gp.exons.insert((300, 400));
    /// assert_eq!(gp.get_exon_count(), 2);
    /// ```
    pub fn get_exon_count(&self) -> usize {
        self.exons.len()
    }

    /// Returns a vector containing the sizes (end coordinates) of all exons.
    ///
    /// Note: This returns the "end" coordinate from the `(start, end)` tuple for each exon.
    /// In genomic contexts, this is often the length or end position, depending on how `exons`
    /// was populated.
    ///
    /// # Returns
    /// A `Vec<u64>` where each element is the end coordinate of an exon.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType}; //
    /// let mut gp = GenePred::new();
    /// gp.exons.insert((10, 20));
    /// gp.exons.insert((30, 40));
    /// let mut sizes = gp.get_exon_sizes();
    /// sizes.sort(); // BTreeSet iteration order is stable, but sorting ensures order for assert
    /// assert_eq!(sizes, vec![20, 40]);
    /// ```
    pub fn get_exon_sizes(&self) -> Vec<u64> {
        self.exons.iter().map(|item| item.1).collect()
    }

    /// Returns a vector containing the start coordinates of all exons.
    ///
    /// # Returns
    /// A `Vec<u64>` where each element is the start coordinate of an exon.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType}; //
    /// let mut gp = GenePred::new();
    /// gp.exons.insert((10, 20));
    /// gp.exons.insert((30, 40));
    /// let mut starts = gp.get_exon_starts();
    /// starts.sort();
    /// assert_eq!(starts, vec![10, 30]);
    /// ```
    pub fn get_exon_starts(&self) -> Vec<u64> {
        self.exons.iter().map(|item| item.0).collect()
    }

    /// Returns a vector containing the relative start coordinates of all exons.
    ///
    /// Each relative start is calculated as `exon_start - gene_pred_start`.
    ///
    /// # Returns
    /// A `Vec<u64>` where each element is the relative start coordinate of an exon.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType};
    /// let mut gp = GenePred::new();
    /// gp.start = 100;
    /// gp.exons.insert((110, 150));
    /// gp.exons.insert((180, 200));
    /// let mut relative_starts = gp.get_exon_starts_relative();
    /// relative_starts.sort();
    /// assert_eq!(relative_starts, vec![10, 80]);
    /// ```
    pub fn get_exon_starts_relative(&self) -> Vec<u64> {
        self.exons
            .iter()
            .map(|item| {
                if self.start > item.0 {
                    panic!(
                        "ERROR: exon start ({}) is less than gene prediction start ({})",
                        item.0, self.start
                    );
                }
                item.0 - self.start
            })
            .collect()
    }

    /// Returns a tuple containing two strings:
    /// 1. A comma-separated string of exon sizes (end coordinates), followed by a comma.
    /// 2. A comma-separated string of relative exon start coordinates, followed by a comma.
    ///
    /// These formats are commonly used in various genomic file specifications.
    ///
    /// # Returns
    /// A tuple `(String, String)` representing `(exon_sizes_string, exon_starts_relative_string)`.
    ///
    /// # Example
    /// ```rust
    /// use std::collections::BTreeSet;
    /// use gxf2bed::{GenePred, Strand, RecordType};
    /// let mut gp = GenePred::new();
    /// gp.start = 100;
    /// gp.exons.insert((110, 150));
    /// gp.exons.insert((180, 200));
    ///
    /// let (exon_sizes, exon_starts_relative) = gp.get_exons_info();
    /// // The order of exons in the BTreeSet is stable (sorted by tuple), so output will be consistent.
    /// assert_eq!(exon_sizes, "150,200,");
    /// assert_eq!(exon_starts_relative, "10,80,");
    /// ```
    pub fn get_exons_info(&self) -> (String, String) {
        let exon_sizes = self
            .get_exon_sizes()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",")
            + ",";

        let exon_starts = self
            .get_exon_starts_relative()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",")
            + ",";

        (exon_sizes, exon_starts)
    }
}

impl<'a> GxfRecord<'a> {
    /// Parses a single line of a GXF (GTF/GFF) file into a `GxfRecord` instance.
    ///
    /// This function expects a tab-separated line and extracts the chromosome, feature,
    /// start, end, strand, frame, and attributes. It also handles parsing the
    /// `attribute` field using the provided `Attribute::parse` method.
    ///
    /// Panics if essential fields (chr, source, feature, start, end, score, strand, frame, attributes)
    /// are missing from the input line.
    ///
    /// # Type Parameters
    /// * `SEP` - A const generic parameter for the separator used within the attribute string.
    ///
    /// # Arguments
    /// * `line` - The input line from the GXF file.
    /// * `attribute` - The name of the attribute key to extract from the attributes column (e.g., "gene_id").
    ///
    /// # Returns
    /// A `Result` containing the parsed `GxfRecord` or a `&'static str` error if the line is empty.
    ///
    /// # Panics
    /// This function will panic if any of the mandatory 9 fields are missing from the input line,
    /// or if the `start` or `end` coordinates cannot be parsed as `u64`.
    ///
    /// # Example
    /// ```rust, ignore
    /// // Assuming Attribute struct and its parse method are available
    /// // and gxf2bed is replaced with the actual crate name.
    /// use gxf2bed::{GxfRecord, Strand, Attribute};
    /// use std::str::FromStr;
    ///
    /// let line = "chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id \"gene1\"; gene_name \"GeneA\";";
    /// let attribute_name = "gene_id".to_string();
    /// let record = GxfRecord::parse::<b';'>(line, &attribute_name).unwrap();
    ///
    /// assert_eq!(record.chr, "chr1");
    /// assert_eq!(record.feature, "gene");
    /// assert_eq!(record.start, 999); // 0-based
    /// assert_eq!(record.end, 5000);
    /// assert_eq!(record.strand, Strand::Forward);
    /// assert_eq!(record.attr.feature(), "gene1"); // Assuming Attribute::parse extracts gene_id as feature
    /// ```
    pub fn parse<const SEP: u8>(line: &'a str, attribute: &String) -> Result<Self, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');

        let (chr, _, feature, start, end, _, strand, frame, attr) = (
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing chromosome from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing source from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing feature from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing start from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing end from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing score from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing strand from GTF records -> {line}"))
                .parse::<Strand>()
                .unwrap_or(Strand::Unknown),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing frame from GTF records -> {line}")),
            fields
                .next()
                .unwrap_or_else(|| panic!("ERROR: Missing attributes from GTF records -> {line}")),
        );

        let attr = Attribute::parse::<SEP>(attr, attribute)
            .map_err(|e| format!("Error parsing attributes: {e}"))
            .unwrap();

        Ok(Self {
            chr: chr.to_string(),
            feature,
            start: start.parse::<u64>().unwrap() - 1,
            end: end.parse().unwrap(),
            strand,
            frame: frame.to_string(),
            attr,
        })
    }

    /// Checks if the `GxfRecord` has a non-empty feature extracted from its attributes.
    ///
    /// This relies on the `Attribute::feature()` method to return the relevant feature string.
    ///
    /// # Returns
    /// `true` if the feature string from attributes is not empty, `false` otherwise.
    ///
    /// # Example
    /// ```rust, ignore
    /// // Assuming Attribute struct and gxf2bed are available
    /// use gxf2bed::{GxfRecord, Strand, Attribute};
    /// use std::str::FromStr;
    ///
    /// let line_with_feature = "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"ABC\";";
    /// let record_with_feature = GxfRecord::parse::<b';'>(line_with_feature, &"gene_id".to_string()).unwrap();
    /// assert!(record_with_feature.has_feature());
    ///
    /// let line_no_feature = "chr1\t.\tgene\t1\t100\t.\t+\t.\tgene_id \"\";";
    /// let record_no_feature = GxfRecord::parse::<b';'>(line_no_feature, &"gene_id".to_string()).unwrap();
    /// assert!(!record_no_feature.has_feature());
    /// ```
    pub fn has_feature(&self) -> bool {
        !self.attr.feature().is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_gtf() {
        let line = "chr1\tunknown\texon\t11869\t12227\t.\t+\t.\tgene_id \"DDX11L1\"; gene_name \"DDX11L1\"; gene_source \"ensembl_havana\";
        gene_biotype \"transcribed_unprocessed_pseudogene\";".to_string();
        let feature = "gene_id".to_string();
        let record = GxfRecord::parse::<b' '>(&line, &feature).unwrap();
        assert_eq!(record.chr, "chr1");
        assert_eq!(record.feature, "exon");
        assert_eq!(record.start, 11868);
        assert_eq!(record.end, 12227);
        assert_eq!(record.strand, Strand::Forward);
        assert_eq!(record.frame, ".");
        assert_eq!(record.attr.feature(), "DDX11L1");
    }

    #[test]
    fn test_record_gff() {
        let line = "chr1\tunknown\texon\t11869\t12227\t.\t+\t.\tID=ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene";
        let feature = "ID".to_string();
        let record = GxfRecord::parse::<b'='>(line, &feature).unwrap();
        assert_eq!(record.chr, "chr1");
        assert_eq!(record.feature, "exon");
        assert_eq!(record.start, 11868);
        assert_eq!(record.end, 12227);
        assert_eq!(record.strand, Strand::Forward);
        assert_eq!(record.frame, ".");
        assert_eq!(record.attr.feature(), "ENSG00000223972");
    }

    #[test]
    fn test_empty_line() {
        let line = "";
        let feature = "ID".to_string();
        let record = GxfRecord::parse::<b' '>(line, &feature);
        assert_eq!(record, Err("Empty line"));
    }

    #[test]
    fn test_empty_strand() {
        let line = "chr1\tunknown\texon\t11869\t12227\t.\t+\t.\tID=ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene";
        let feature = "ID".to_string();
        let record = GxfRecord::parse::<b'='>(line, &feature).unwrap();
        assert_eq!(record.strand, Strand::Forward);
    }

    #[test]
    fn test_gene_pred() {
        let mut gene_pred = GenePred::new();
        let query = GenePred {
            chr: "chr1".to_string(),
            start: 11868,
            end: 12227,
            strand: Strand::Forward,
            exons: vec![(11868, 50), (12200, 100)].into_iter().collect(),
            record_type: RecordType::Parent,
            blocks: None,
        };

        gene_pred.merge(query);

        assert_eq!(gene_pred.chr, "chr1");
        assert_eq!(gene_pred.start, 11868);
        assert_eq!(gene_pred.end, 12227);
        assert_eq!(gene_pred.strand, Strand::Forward);
        assert_eq!(gene_pred.get_exon_count(), 2);
        assert_eq!(gene_pred.get_exon_sizes(), vec![50, 100]);
        assert_eq!(gene_pred.get_exon_starts(), vec![11868, 12200]);
        assert_eq!(gene_pred.get_exon_starts_relative(), vec![0, 332]);
        assert_eq!(
            gene_pred.get_exons_info(),
            ("50,100,".to_string(), "0,332,".to_string())
        );
    }
}
