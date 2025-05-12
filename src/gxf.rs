mod attr;
pub use attr::*;

use std::collections::BTreeSet;

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

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Strand::Forward => write!(f, "+"),
            Strand::Reverse => write!(f, "-"),
            Strand::Unknown => write!(f, "."),
        }
    }
}

#[derive(Debug, PartialEq)]
pub enum RecordType {
    Parent,
    Child,
    Unknown,
}

#[derive(Debug, PartialEq)]
pub struct GenePred {
    pub chr: String,
    pub start: u64,
    pub end: u64,
    pub strand: Strand,
    pub exons: BTreeSet<(u64, u64)>,
    pub record_type: RecordType,
    pub cds_start: Option<u64>,
    pub cds_end: Option<u64>,
}

impl GenePred {
    pub fn new() -> Self {
        Self {
            chr: String::new(),
            start: 0,
            end: 0,
            strand: Strand::Unknown,
            exons: BTreeSet::new(),
            record_type: RecordType::Unknown,
            cds_start: None,
            cds_end: None,
        }
    }

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

        if let Some(qs) = query.cds_start {
            self.cds_start = Some(self.cds_start.map_or(qs, |cs| cs.min(qs)));
        }
        if let Some(qe) = query.cds_end {
            self.cds_end = Some(self.cds_end.map_or(qe, |ce| ce.max(qe)));
        }
    }

    pub fn get_exon_count(&self) -> usize {
        self.exons.len()
    }

    pub fn get_exon_sizes(&self) -> Vec<u64> {
        self.exons.iter().map(|item| item.1).collect()
    }

    pub fn get_exon_starts(&self) -> Vec<u64> {
        self.exons.iter().map(|item| item.0).collect()
    }

    pub fn get_exon_starts_relative(&self) -> Vec<u64> {
        self.exons
            .iter()
            .map(|item| {
                if self.start > item.0 {
                    dbg!(&self);
                }
                item.0 - self.start
            })
            .collect()
    }

    pub fn get_cds_start(&self) -> u64 {
        self.cds_start.unwrap_or(self.start)
    }

    pub fn get_cds_end(&self) -> u64 {
        self.cds_end.unwrap_or(self.end)
    }

    pub fn get_cds(&self) -> (u64, u64) {
        (self.get_cds_start(), self.get_cds_end())
    }

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
    pub fn parse<const SEP: u8>(line: &'a str, attribute: &String) -> Result<Self, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');

        let (chr, _, feature, start, end, _, strand, frame, attr) = (
            fields.next().ok_or("Missing chrom")?,
            fields.next().ok_or("Missing source")?,
            fields.next().ok_or("Missing feature")?,
            fields.next().ok_or("Missing start")?,
            fields.next().ok_or("Missing end")?,
            fields.next().ok_or("Missing score")?,
            fields.next().ok_or("Missing strand")?,
            fields.next().ok_or("Missing frame")?,
            fields.next().ok_or("Missing attributes")?,
        );

        let strand = match strand.chars().next().expect("ERROR: Strand is empty") {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => Strand::Unknown,
        };

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
            cds_start: Some(11868),
            cds_end: Some(12300),
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
        assert_eq!(gene_pred.get_cds_start(), 11868);
        assert_eq!(gene_pred.get_cds_end(), 12300);
        assert_eq!(gene_pred.get_cds(), (11868, 12300));
        assert_eq!(
            gene_pred.get_exons_info(),
            ("50,100,".to_string(), "0,332,".to_string())
        );
    }
}
