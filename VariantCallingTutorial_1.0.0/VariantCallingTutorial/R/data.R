.load_vcf_path <- function(x) {
  TabixFile(system.file("extdata", x, package="VariantCallingTutorial"))
}
