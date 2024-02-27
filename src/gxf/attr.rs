// use std::collections::HashMap;
use hashbrown::HashMap;
use std::fmt::Debug;

use thiserror::Error;

#[derive(Debug, PartialEq)]
pub struct Attribute {
    feature: String,
}

impl Attribute {
    pub fn parse(line: &String, feature: &String) -> Result<Self, ParseError> {
        let mut attributes: HashMap<String, String> = HashMap::new();
        let bytes = line.trim_end().as_bytes().iter().enumerate();

        let mut start = 0;

        for (mut i, byte) in bytes {
            if *byte == b';' || i == line.len() - 1 {
                if i == line.len() - 1 && *byte != b';' {
                    i += 1;
                };
                let word = &line[start..i];
                if !word.is_empty() {
                    let (key, value) = get_pair(word).ok_or(ParseError::Invalid)?;
                    attributes.insert(key, value);
                }
                start = i + 1;
            }
        }

        Ok(Attribute {
            feature: attributes
                .get(feature)
                .unwrap_or(&"".to_string())
                .to_string(),
        })
    }

    pub fn feature(&self) -> &str {
        &self.feature
    }
}

fn get_pair(line: &str) -> Option<(String, String)> {
    let line = line.trim();
    let mut bytes = line.as_bytes().iter();
    let i = bytes
        .position(|b| *b == b' ' || *b == b'=')
        .ok_or(ParseError::Invalid)
        .map_err(|e| {
            eprintln!("{:?}", e);
            e
        })
        .unwrap();

    let key = &line[..i];
    let value = &line[i + 1..line.len()].trim_matches('"').trim();

    Some((key.to_string(), value.to_string()))
}

#[derive(Error, Debug, PartialEq)]
pub enum ParseError {
    #[error("Empty line")]
    Empty,
    #[error("Invalid GTF line")]
    Invalid,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_pair_gtf() {
        let line = " gene_name \"DDX11L1\"";
        let (key, value) = get_pair(line).unwrap();
        assert_eq!(key, "gene_name");
        assert_eq!(value, "DDX11L1");
    }

    #[test]
    fn test_parse_gtf() {
        let line = "gene_id \"ENSG00000223972\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene OTTHUMG00000000961.1;";
        let feature = "gene_id".to_string();
        let attr = Attribute::parse(&line.to_string(), &feature).unwrap();
        assert_eq!(attr.feature, "ENSG00000223972");
    }

    #[test]
    fn test_get_pair_gff() {
        let line = "ID=ENSG00000223972";
        let (key, value) = get_pair(line).unwrap();
        assert_eq!(key, "ID");
        assert_eq!(value, "ENSG00000223972");
    }

    #[test]
    fn test_parse_gff() {
        let line = "ID=ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene";
        let feature = "ID".to_string();
        let attr = Attribute::parse(&line.to_string(), &feature).unwrap();
        assert_eq!(attr.feature, "ENSG00000223972");
    }
}
