library(httr); library(jsonlite); library(data.table)
ot_endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

## ---- 1) 先拿 tractability（这个字段确认可用，先保证有产出）----
q_tract <- '
{
  PTK2:  target(ensemblId:"ENSG00000169398"){ approvedSymbol tractability{ label modality value } }
  ROCK1: target(ensemblId:"ENSG00000067900"){ approvedSymbol tractability{ label modality value } }
  DLK1:  target(ensemblId:"ENSG00000185559"){ approvedSymbol tractability{ label modality value } }
  SIK2:  target(ensemblId:"ENSG00000170145"){ approvedSymbol tractability{ label modality value } }
}'
r1 <- POST(ot_endpoint, body=list(query=q_tract), encode="json")
cat("tractability HTTP:", status_code(r1), "\n")
d1 <- content(r1, as="parsed", simplifyVector=FALSE)
if(!is.null(d1$errors)) { cat("errors:\n"); str(d1$errors, max.level=2) }
for(g in c("PTK2","ROCK1","DLK1","SIK2")){
  t <- d1$data[[g]]; if(is.null(t)) next
  tr <- rbindlist(lapply(t$tractability, function(x)
        data.table(modality=x$modality, label=x$label, value=isTRUE(x$value))), fill=TRUE)
  cat(sprintf("\n%s — tractability (value=TRUE):\n", g))
  print(tr[value==TRUE, .(modality, label)])
}

## ---- 2) 内省 Target 类型的所有字段，找“药物”字段的真名 ----
q_intro <- '{ __type(name:"Target"){ fields{ name } } }'
r2 <- POST(ot_endpoint, body=list(query=q_intro), encode="json")
d2 <- content(r2, as="parsed", simplifyVector=FALSE)
fields <- sapply(d2$data[["__type"]]$fields, function(x) x$name)
cat("\n=== Target 类型全部字段（共", length(fields), "个）===\n")
print(fields)
cat("\n>>> 药物相关字段（drug/chembl 等）:\n")
print(grep("drug|Drug|chembl|Chembl|medicine|Medicine", fields, value=TRUE))

library(httr); library(jsonlite); library(data.table)
ot_endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

## ---- 1) 内省 drugAndClinicalCandidates 的返回结构，拿到正确子字段名 ----
q_intro2 <- '
{
  __type(name:"Target"){
    fields(includeDeprecated:true){
      name
      args{ name }
      type{ name kind ofType{ name kind ofType{ name kind } } }
    }
  }
}'
r <- POST(ot_endpoint, body=list(query=q_intro2), encode="json")
d <- content(r, as="parsed", simplifyVector=FALSE)
fld <- Filter(function(x) x$name=="drugAndClinicalCandidates", d$data[["__type"]]$fields)[[1]]
cat("drugAndClinicalCandidates 的参数:\n"); print(sapply(fld$args, function(a) a$name))
## 找到它返回的对象类型名，再内省该类型的字段
tn <- fld$type$name; if(is.null(tn)) tn <- fld$type$ofType$name; if(is.null(tn)) tn <- fld$type$ofType$ofType$name
cat("返回对象类型:", tn, "\n")
q_intro3 <- sprintf('{ __type(name:"%s"){ fields{ name type{ name kind ofType{ name } } } } }', tn)
d3 <- content(POST(ot_endpoint, body=list(query=q_intro3), encode="json"), as="parsed", simplifyVector=FALSE)
cat("该类型的字段:\n"); print(sapply(d3$data[["__type"]]$fields, function(x) x$name))

library(httr); library(jsonlite); library(data.table)
ot_endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"

## 标准内省：把 clinicalTargets 的字段和完整 ofType 链一次拿全
q <- '
{
  __type(name:"clinicalTargets"){
    fields{
      name
      type{ kind name ofType{ kind name ofType{ kind name ofType{ kind name ofType{ kind name } } } } }
    }
  }
}'
d <- content(POST(ot_endpoint, body=list(query=q), encode="json"), as="parsed", simplifyVector=FALSE)
rowsfld <- Filter(function(x) x$name=="rows", d$data[["__type"]]$fields)[[1]]

## 递归钻 ofType，直到拿到底层的具名类型（OBJECT）
dig <- function(t){ while(!is.null(t) && is.null(t$name)) t <- t$ofType; if(is.null(t)) return(NA); t$name }
rt <- dig(rowsfld$type)
cat("rows 元素类型:", rt, "\n")

## 列出该类型所有字段
if(!is.na(rt)){
  q2 <- sprintf('{ __type(name:"%s"){ fields{ name type{ kind name ofType{ kind name ofType{ kind name } } } } } }', rt)
  d2 <- content(POST(ot_endpoint, body=list(query=q2), encode="json"), as="parsed", simplifyVector=FALSE)
  cat("\n=== rows 每条记录的字段 ===\n")
  for(x in d2$data[["__type"]]$fields){
    tn <- dig(x$type)
    cat(sprintf("  %-24s : %s\n", x$name, ifelse(is.na(tn),"",tn)))
  }
}

library(httr); library(jsonlite); library(data.table)
ot_endpoint <- "https://api.platform.opentargets.org/api/v4/graphql"
targets <- c(PTK2="ENSG00000169398", ROCK1="ENSG00000067900",
             DLK1="ENSG00000185559", SIK2="ENSG00000170145")

## ---- 1) 取 drugAndClinicalCandidates（真名：maxClinicalStage / drug{...} / diseases）----
build_q <- function(alias, ensg) sprintf('
  %s: target(ensemblId:"%s"){
    approvedSymbol
    drugAndClinicalCandidates{
      count
      rows{
        maxClinicalStage
        drug{ id name isApproved drugType }
        diseases{ name }
      }
    }
  }', alias, ensg)
q <- paste0("{\n", paste(mapply(build_q, names(targets), targets), collapse="\n"), "\n}")
d <- content(POST(ot_endpoint, body=list(query=q), encode="json"), as="parsed", simplifyVector=FALSE)
if(!is.null(d$errors)){ cat("errors:\n"); str(d$errors, max.level=2) }

## 解析每个靶：药物数 + 各药(阶段/批准/适应症数)
ot_drug_list <- list()
for(g in names(targets)){
  t <- d$data[[g]]; dc <- t$drugAndClinicalCandidates
  cat(sprintf("\n===== %s ===== 候选药条目: %s\n", g, if(is.null(dc)) 0 else dc$count))
  if(is.null(dc) || dc$count==0){ next }
  rows <- rbindlist(lapply(dc$rows, function(r) data.table(
    Gene       = g,
    drug       = if(is.null(r$drug$name)) NA_character_ else r$drug$name,
    drugType   = if(is.null(r$drug$drugType)) NA_character_ else r$drug$drugType,
    approved   = isTRUE(r$drug$isApproved),
    maxStage   = if(is.null(r$maxClinicalStage)) NA_character_ else r$maxClinicalStage,
    n_disease  = length(r$diseases))), fill=TRUE)
  ot_drug_list[[g]] <- rows
  # 按药去重展示（阶段最高在前）
  stage_ord <- c("Approved"=5,"Phase IV"=4,"Phase III"=3,"Phase II"=2,"Phase I"=1)
  rows[, stg_rank := stage_ord[maxStage]]
  print(unique(rows[order(-stg_rank)], by="drug")[, .(drug, drugType, approved, maxStage)][1:min(8,.N)])
}
ot_drugs <- rbindlist(ot_drug_list, fill=TRUE)
fwrite(ot_drugs, "E:/pulp/druggability_opentargets_drugs.csv")

## ---- 2) 汇总每靶：最高阶段 + 已批准药数 + 已批准药名 ----
ot_summary <- ot_drugs[, .(
  OT_n_candidates = uniqueN(drug),
  OT_max_stage    = maxStage[which.max(c("Approved"=5,"Phase IV"=4,"Phase III"=3,"Phase II"=2,"Phase I"=1)[maxStage])],
  OT_approved_drugs = paste(unique(drug[approved==TRUE]), collapse="; ")
), by=Gene]
cat("\n=== OT 每靶汇总 ===\n"); print(ot_summary)


rlibrary(data.table)
## 用 OT 实际返回的阶段字符串重建映射
stage_rank <- c("APPROVAL"=5, "PHASE_4"=4, "PHASE_3"=3, "PHASE_2_3"=2.5,
                "PHASE_2"=2, "PHASE_1_2"=1.5, "PHASE_1"=1, "EARLY_PHASE_1"=0.5)
ot_drugs <- fread("E:/pulp/druggability_opentargets_drugs.csv")

ot_summary <- ot_drugs[, {
  u <- unique(.SD, by="drug")
  .(OT_n_drugs   = nrow(u),
    OT_max_stage = u$drugMaxStage[which.max(stage_rank[u$drugMaxStage])],
    OT_approved  = paste(unique(u$drug[u$drugMaxStage=="APPROVAL"]), collapse="; "))
}, by=Gene]
cat("=== OT 每靶汇总（修正阶段映射）===\n"); print(ot_summary)
跑完 OT_max_stage 应显示 ROCK1=APPROVAL、PTK2=PHASE_3,OT_approved 里 ROCK1 应列出 FASUDIL; NETARSUDIL; BELUMOSUDIL; RIPASUDIL。