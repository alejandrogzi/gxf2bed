mod attr;
pub use attr::*;

#[derive(Debug, PartialEq)]
pub struct GxfRecord<'a> {
    pub chr: String,
    pub feat: String,
    pub start: u32,
    pub end: u32,
    pub strand: char,
    pub frame: String,
    pub attr: Attribute<'a>,
}

impl<'a> GxfRecord<'a> {
    pub fn parse<const SEP: u8>(line: &'a str, feature: &String) -> Result<Self, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let mut fields = line.split('\t');

        let (chr, _, feat, start, end, _, strand, frame, attr) = (
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

        let attr = Attribute::parse::<SEP>(attr, feature)
            .map_err(|e| format!("Error parsing attributes: {e}"))
            .unwrap();

        Ok(Self {
            chr: chr.to_string(),
            feat: feat.to_string(),
            start: start.parse::<u32>().unwrap() - 1,
            end: end.parse().unwrap(),
            strand: strand.chars().next().unwrap(),
            frame: frame.to_string(),
            attr,
        })
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
        assert_eq!(record.feat, "exon");
        assert_eq!(record.start, 11868);
        assert_eq!(record.end, 12227);
        assert_eq!(record.strand, '+');
        assert_eq!(record.frame, ".");
        assert_eq!(record.attr.feature(), "DDX11L1");
    }

    #[test]
    fn test_record_gff() {
        let line = "chr1\tunknown\texon\t11869\t12227\t.\t+\t.\tID=ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene";
        let feature = "ID".to_string();
        let record = GxfRecord::parse::<b'='>(line, &feature).unwrap();
        assert_eq!(record.chr, "chr1");
        assert_eq!(record.feat, "exon");
        assert_eq!(record.start, 11868);
        assert_eq!(record.end, 12227);
        assert_eq!(record.strand, '+');
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
}
