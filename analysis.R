# Run installation only once on your computer:
install.packages(c("rentrez", "ggplot2"))
library(rentrez)
library(ggplot2)
# CHANGETHESEFOURVALUEStoyourassignment (example: ACTA1):
gene_symbol <- "MYH2"
transcript_accession <- "NM_017534.6"
clinvar_id <- "625152"
variant_protein <- "p.Arg793Ter" # Use "" for splice variants
for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}

hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)
if (length(hits$ids) == 0) stop("Gene not found")
# Verify the exact gene symbol, not just a similar family member.
gene_id <- NA_character_
for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id; gene_info <- rec; break
  }
}
if (is.na(gene_id)) stop("Exact human gene not found")
print(gene_info)
writeLines(c(
  paste("Gene symbol:", gene_info$name),
  paste("NCBI Gene ID:", gene_id),
  paste("Official name:", gene_info$description),
  paste("Chromosome:", gene_info$chromosome)
), "results/gene_summary.txt")

# Display and verify the retrieved gene information
print(gene_info$name)
print(gene_info$description)
print(gene_info$chromosome)
print(gene_id)

# Check the saved summary file
readLines("results/gene_summary.txt")

# Add the normal role of MYH2 in muscle
gene_role <- paste(
  "MYH2 encodes myosin heavy chain 2, a motor protein",
  "found in fast-twitch skeletal muscle fibers.",
  "It interacts with actin and uses ATP to generate force",
  "during muscle contraction."
)

cat(gene_role, "\n")

write(
  gene_role,
  file = "results/gene_summary.txt",
  append = TRUE
)


# Step 3: Retrieve the matching protein sequence

gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0 || fragment == "") {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub('^/protein_id="|"$', '', fragment)

protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA")
}

fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

cat("Matched protein:", protein_accession, "\n")

writeLines(
  paste("NP accession:", protein_accession),
  "results/protein_accession.txt"
)


# Step 4: Analyze protein length and amino-acid composition

fasta_lines <- readLines(fasta_file)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]
protein_length <- length(aa)

standard <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
counts <- table(factor(aa, levels = standard))

results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

cat("Protein length:", protein_length, "amino acids\n")

print(head(results[order(-results$Count), ], 3))

fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = paste(gene_symbol, "amino-acid composition"),
    x = "Amino acid",
    y = "Number of residues"
  )

print(fig1)

ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)


# Step 6: Improved MYH2 variant-position plot

position <- 793

fig2 <- ggplot() +
  # Full reference protein
  geom_segment(
    aes(x = 1, xend = protein_length, y = 0, yend = 0),
    linewidth = 3,
    color = "grey65"
  ) +
  
  # Variant position
  geom_point(
    aes(x = position, y = 0),
    color = "red",
    size = 4
  ) +
  
  # Label above the variant
  annotate(
    "label",
    x = position,
    y = 0.22,
    label = "p.Arg793Ter\nPosition 793",
    color = "red",
    size = 3.5
  ) +
  
  # Start and end labels
  annotate("text", x = 1, y = -0.18, label = "1", size = 4) +
  annotate(
    "text",
    x = protein_length,
    y = -0.18,
    label = "1,941 aa",
    hjust = 1,
    size = 4
  ) +
  
  scale_x_continuous(
    limits = c(1, protein_length),
    breaks = c(1, 400, 793, 1200, 1600, 1941),
    labels = c("1", "400", "793", "1200", "1600", "1941")
  ) +
  
  coord_cartesian(ylim = c(-0.35, 0.45), clip = "off") +
  
  labs(
    title = "MYH2 Variant Position on the Reference Protein",
    subtitle = "Reference protein NP_060004.3 | Total length: 1,941 amino acids",
    x = "Amino-acid position",
    y = NULL,
    caption = paste(
      "The red marker identifies the stop-gain variant at residue 793.",
      "It does not indicate a confirmed truncated-protein length."
    )
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8),
    plot.margin = margin(15, 20, 15, 20)
  )

print(fig2)

ggsave(
  "figures/02_variant_position.png",
  fig2,
  width = 10,
  height = 4,
  dpi = 300
)



# Step 8: Improved MYH2 functional-domain map

fig3 <- ggplot() +
  # Full reference protein
  geom_segment(
    aes(x = 1, xend = protein_length, y = 0, yend = 0),
    linewidth = 2,
    color = "grey75"
  ) +
  
  # Annotated domains
  geom_rect(
    data = domain_data,
    aes(
      xmin = start, xmax = end,
      ymin = -0.10, ymax = 0.10,
      fill = domain
    ),
    color = "black",
    linewidth = 0.3
  ) +
  
  # Variant marker
  geom_vline(
    xintercept = 793,
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "point", x = 793, y = 0,
    color = "red", size = 3
  ) +
  annotate(
    "label",
    x = 793, y = 0.28,
    label = "p.Arg793Ter\nPosition 793",
    color = "red",
    size = 3.5
  ) +
  
  # Avoid crowded tick labels
  scale_x_continuous(
    breaks = c(1, 100, 400, 600, 1000, 1200, 1600, 1941),
    labels = c("1", "100", "400", "600",
               "1000", "1200", "1600", "1941"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  scale_fill_brewer(palette = "Set2") +
  
  coord_cartesian(
    xlim = c(1, protein_length),
    ylim = c(-0.20, 0.45),
    clip = "off"
  ) +
  
  labs(
    title = "MYH2 Functional-Domain Map",
    subtitle = paste(
      "Reference:", protein_accession,
      "| Length:", protein_length, "amino acids"
    ),
    x = "Amino-acid position",
    y = NULL,
    fill = "Annotated domain",
    caption = paste(
      "Domain coordinates based on NCBI MYH2 conserved-domain annotations.",
      "Red marker: p.Arg793Ter at residue 793."
    )
  ) +
  
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8),
    plot.margin = margin(12, 15, 12, 15)
  )

print(fig3)

ggsave(
  "figures/03_protein_structure.png",
  fig3,
  width = 10,
  height = 4,
  dpi = 300
)
