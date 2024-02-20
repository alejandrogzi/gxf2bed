use std::collections::HashMap;
use std::fmt::Debug;

use thiserror::Error;

#[derive(Debug)]
pub struct Attribute {
    // gene_id: String,
    feature: String,
}

impl Attribute {
    pub fn parse(line: &String, feature: &String) -> Result<Self, ParseError> {
        let mut attributes: HashMap<String, String> = HashMap::new();
        let bytes = line.as_bytes().iter().enumerate();

        let mut start = 0;

        for (i, byte) in bytes {
            if *byte == b';' {
                let word = &line[start..i];
                if !word.is_empty() {
                    let (key, value) = get_pair(word).ok_or(ParseError::Invalid)?;
                    attributes.insert(key, value);
                }
                start = i + 1;
            }
        }

        // let gene_id = attributes.get("gene_id").ok_or(ParseError::Invalid);
        Ok(Attribute {
            // gene_id: gene_id?.to_string(),
            feature: attributes
                .get(feature)
                .unwrap_or(&"".to_string())
                .to_string(),
        })
    }
    //
    // pub fn gene_id(&self) -> &str {
    //     &self.gene_id
    // }

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
