.onLoad <- function(libname, pkgName) {
  ns <- asNamespace(pkgName)

  delayedAssign("NA12878_pg.chr20.vcf.bgz",
                .load_vcf_path("NA12878_pg.chr20.vcf.bgz"),
                eval.env = ns, assign.env = ns)
  namespaceExport(ns, ls(ns)) 
}
