write_chr_pos_metric <- function(
  chr, pos, value, file, value_name
) {
  data.table::fwrite(
    data.table::data.table(
      Chromosome = chr,
      Position   = pos,
      value      = value
    ),
    file = file,
    sep = "\t",
    col.names = c("Chromosome", "Position", value_name)
  )
}
