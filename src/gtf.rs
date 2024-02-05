mod attr;
pub use attr::*;

#[derive(Debug)]
pub struct GtfRecord {
    pub chr: String,
    pub feat: String,
    pub start: u32,
    pub end: u32,
    pub strand: char,
    pub frame: String,
    pub attr: Attribute,
}

impl GtfRecord {
    pub fn parse(line: &str) -> Result<Self, &'static str> {
        if line.is_empty() {
            return Err("Empty line");
        }

        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 9 {
            return Err("Line has fewer fields than expected".into());
        }

        let attr =
            Attribute::parse(&fields[8].to_string()).map_err(|_| "Error parsing attribute")?;

        Ok(Self {
            chr: String::from(fields[0]),
            feat: String::from(fields[2]),
            start: fields[3].parse::<u32>().unwrap() - 1,
            end: fields[4].parse().unwrap(),
            strand: fields[6].parse().unwrap(),
            frame: String::from(fields[7]),
            attr: attr,
        })
    }
}
