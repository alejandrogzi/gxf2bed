//! The fastest GTF/GFF-to-BED converter chilling around
//! Alejandro Gonzales-Irribarren, 2025

use std::fmt::Debug;
use thiserror::Error;

/// A macro to extract a field's value from a byte slice.
///
/// This macro provides two rules for extracting key-value pairs from a byte slice,
/// typically representing attributes in a GXF file. It supports both a variable separator
/// and a literal byte separator for the key-value pairs.
///
/// The macro attempts to strip the `field_name` prefix, then the separator `sep`,
/// and finally trims surrounding double quotes from the value. The extracted value
/// is converted to a `&str` using `unsafe { std::str::from_utf8_unchecked(without_eq) }`,
/// assuming valid UTF-8 input after stripping.
///
/// # Usage
/// ```rust, ignore
/// let bytes_data = b"gene_id \"gene1\";";
/// let mut gene_id_val = None;
/// let separator = b' ';
/// extract_field!(
///     bytes_data split by separator to
///     b"gene_id" => &mut gene_id_val;
/// );
/// assert_eq!(gene_id_val, Some("gene1"));
///
/// let bytes_data_gff = b"ID=gene1;";
/// let mut id_val = None;
/// extract_field!(
///     bytes_data_gff split by b'=' to
///     b"ID" => &mut id_val;
/// );
/// assert_eq!(id_val, Some("gene1"));
/// ```
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

/// Splits a byte slice by a specified byte separator and trims a leading byte from each resulting sub-slice.
///
/// This function is particularly useful for parsing attribute strings in GXF files
/// where individual attributes are separated by a delimiter (e.g., `;`) and might
/// have leading whitespace (` `).
///
/// # Type Parameters
/// * `BY` - The byte to split the slice by (e.g., `b';'`).
/// * `TRIM` - The leading byte to trim from each resulting sub-slice (e.g., `b' '`).
///
/// # Arguments
/// * `bytes` - The input byte slice to split and trim.
///
/// # Returns
/// An iterator over byte slices, where each slice is a segment from the original input,
/// split by `BY` and with leading `TRIM` bytes removed.
///
/// # Example
/// ```rust
/// use gxf2bed::split_and_trim_bytes; //
///
/// let data = b" attr1=val1;  attr2=val2; attr3=val3";
/// let result: Vec<&[u8]> = split_and_trim_bytes::<b';', b' '>(data).collect();
///
/// assert_eq!(result, vec![
///     b"attr1=val1",
///     b"attr2=val2",
///     b"attr3=val3"
/// ]);
///
/// let empty_data = b"";
/// let result_empty: Vec<&[u8]> = split_and_trim_bytes::<b';', b' '>(empty_data).collect();
/// assert!(result_empty.is_empty());
/// ```
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

/// Represents the attributes parsed from a GXF record.
///
/// This struct currently focuses on extracting a single "feature" attribute
/// (e.g., `gene_id` or `transcript_id`) from the attribute string.
///
/// # Fields
/// * `feature` - The value of the extracted feature attribute (e.g., gene ID).
#[derive(Debug, PartialEq)]
pub struct Attribute<'a> {
    feature: &'a str,
}

impl<'a> Attribute<'a> {
    /// Parses an attribute string from a GXF record into an `Attribute` struct.
    ///
    /// This method uses the `split_and_trim_bytes` function and the `extract_field!` macro
    /// to efficiently parse key-value pairs from the attribute string. It specifically
    /// looks for the attribute identified by the `feature` parameter (e.g., "gene_id").
    ///
    /// # Type Parameters
    /// * `SEP` - The byte separator used between the key and value within each attribute (e.g., `b' '` for GTF, `b'='` for GFF).
    ///
    /// # Arguments
    /// * `line` - The attribute string portion of a GXF record.
    /// * `feature` - The name of the attribute key to extract (e.g., "gene_id", "transcript_id").
    ///
    /// # Returns
    /// A `Result` containing the parsed `Attribute` or a `ParseError` if the line is empty.
    /// Note: If the `feature` attribute is not found, `feature()` will return an empty string.
    ///
    /// # Example
    /// ```rust
    /// use gxf2bed::{Attribute, ParseError}; //
    ///
    /// // Example for GTF-like attributes (space separated key-value)
    /// let gtf_attr_line = "gene_id \"gene1\"; transcript_id \"tx1\";";
    /// let attr = Attribute::parse::<b' '>(gtf_attr_line, &"gene_id".to_string()).unwrap();
    /// assert_eq!(attr.feature(), "gene1");
    ///
    /// // Example for GFF-like attributes (equal sign separated key-value)
    /// let gff_attr_line = "ID=gene1;Name=GeneA;";
    /// let attr_gff = Attribute::parse::<b'='>(gff_attr_line, &"ID".to_string()).unwrap();
    /// assert_eq!(attr_gff.feature(), "gene1");
    ///
    /// // Example with missing feature
    /// let missing_feature_line = "other_attr \"val\";";
    /// let attr_missing = Attribute::parse::<b' '>(missing_feature_line, &"gene_id".to_string()).unwrap();
    /// assert_eq!(attr_missing.feature(), "");
    ///
    /// // Example with empty line
    /// let empty_line = "";
    /// let error = Attribute::parse::<b' '>(empty_line, &"gene_id".to_string()).unwrap_err();
    /// assert_eq!(error, ParseError::Empty);
    /// ```
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

    /// Returns the extracted feature string.
    ///
    /// This is an `inline(always)` function, suggesting it's very lightweight and
    /// will likely be inlined by the compiler for performance.
    ///
    /// # Returns
    /// A string slice (`&'a str`) representing the feature value. Returns an empty
    /// string if the feature was not found during parsing.
    ///
    /// # Example
    /// ```rust
    /// use gxf2bed::Attribute; //
    ///
    /// let attr = Attribute { feature: "gene_A" };
    /// assert_eq!(attr.feature(), "gene_A");
    /// ```
    #[inline(always)]
    pub fn feature(&self) -> &'a str {
        self.feature
    }
}

/// Represents errors that can occur during parsing of GXF attributes.
///
/// This enum uses `thiserror` attributes for ergonomic error handling and display.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum ParseError {
    /// Indicates that the input line or attribute string was empty.
    #[error("Empty line, cannot parse attributes")]
    Empty,

    /// Indicates that a specific attribute (e.g., `gene_id`) was missing from the record.
    /// Contains the string of the record where the attribute was missing.
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
