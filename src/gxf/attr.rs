// use hashbrown::HashMap;
use std::fmt::Debug;
use thiserror::Error;

macro_rules! extract_field {
    ($bytes:ident split by $sep:ident to $( $field_name:expr => $output_field:expr; )+) => {
        $(
            if let Some(without_key) = $bytes.strip_prefix($field_name) {
                if let Some(without_eq) = without_key.strip_prefix(&[$sep]) {
                    let value = unsafe { std::str::from_utf8_unchecked(without_eq) };
                    *$output_field = Some(value.trim_matches(|c| c == '"'));
                }
            }
        )+
    };
    ($bytes:ident split by $sep:literal to $( $field_name:literal => $output_field:expr; )+) => {
        $(
            if let Some(without_key) = $bytes.strip_prefix($field_name) {
                if let Some(without_eq) = without_key.strip_prefix(&[$sep]) {
                    let value = unsafe { std::str::from_utf8_unchecked(without_eq) };
                    *$output_field = Some(value.trim_matches(|c| c == '"'));
                }
            }
        )+
    };
}

#[inline(always)]
fn split_and_trim_bytes<const BY: u8, const TRIM: u8>(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    bytes.split(|b| *b == BY).map(|b| {
        let mut idx = 0;
        while idx < b.len() && b[idx] == TRIM {
            idx += 1;
        }
        &b[idx..]
    })
}

#[derive(Debug, PartialEq)]
pub struct Attribute<'a> {
    feature: &'a str,
}

impl<'a> Attribute<'a> {
    pub fn parse<const SEP: u8>(
        line: &'a str,
        feature: &String,
    ) -> Result<Attribute<'a>, ParseError> {
        if !line.is_empty() {
            let field_bytes = split_and_trim_bytes::<b';', b' '>(line.trim_end().as_bytes());

            let mut feat = None;

            for field in field_bytes {
                extract_field!(
                    field split by SEP to
                    feature.as_bytes() => &mut (feat);
                )
            }

            Ok(Attribute {
                feature: feat.unwrap_or(""),
            })
        } else {
            Err(ParseError::Empty)
        }
    }

    #[inline(always)]
    pub fn feature(&self) -> &'a str {
        self.feature
    }
}

#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseError {
    // Empty line
    #[error("Empty line, cannot parse attributes")]
    Empty,

    // Missing gene_id attribute
    #[error("Missing attribute in: {0}")]
    MissingFeature(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gtf() {
        let line = "gene_id \"ENSG00000223972\"; gene_type \"transcribed_unprocessed_pseudogene\"; gene_name \"DDX11L1\"; level 2; havana_gene OTTHUMG00000000961.1;";
        let feature = "gene_id".to_string();
        let attr = Attribute::parse::<b' '>(&line, &feature).unwrap();
        assert_eq!(attr.feature, "ENSG00000223972");
    }

    #[test]
    fn test_parse_gff() {
        let line = "ID=ENSG00000223972;Name=DDX11L1;biotype=transcribed_unprocessed_pseudogene";
        let feature = "ID".to_string();
        let attr = Attribute::parse::<b'='>(&line, &feature).unwrap();
        assert_eq!(attr.feature, "ENSG00000223972");
    }
}
