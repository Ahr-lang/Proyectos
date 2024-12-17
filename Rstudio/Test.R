#A00838626 Andres Huerta Robinson
# A00840653 German Mounier Zalenski
#Victor Valero a01383804
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")
if (!requireNamespace("phangorn", quietly = TRUE)) install.packages("phangorn")
if (!requireNamespace("phytools", quietly = TRUE)) install.packages("phytools")
if (!requireNamespace("geiger", quietly = TRUE)) install.packages("geiger")
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings")
if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr")
if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggtree", quietly = TRUE)) BiocManager::install("ggtree")
if (!requireNamespace("DECIPHER", quietly = TRUE)) BiocManager::install("DECIPHER")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggmsa", quietly = TRUE)) install.packages("ggmsa")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("stringdist", quietly = TRUE)) install.packages("stringdist")
if (!requireNamespace("proxy", quietly = TRUE)) install.packages("proxy")

library(RColorBrewer)
library(stringdist)  # for calculating string distances
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(rentrez)
library(Biostrings)
library(seqinr)
library(adegenet)
library(tidyr)  # Añadir esta línea para cargar tidyr
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggplot2) #forma parte de las paqueterias que comprenden el "nucleo principal" de tidyverse
library(ggmsa)
library(dplyr)
library(adegenet)
#Para esta seccion del proyecto, analizare la proteina N, osea la proteina fosfoproteína nucleocápside.
#Esta proteina  del virus desempeña roles esenciales en la replicación del ARN, el ensamblaje de nuevas partículas virales, y manipula las defensas celulares para favorecer la infección(Li et al., 2020).
#La escogi debido a que me intereso su papel como la estructura de su ARN, pero tambien el que de ella dependa pasar la replicacion, osea no la hace como tal pero su rol es mas como el de guardar y mandar los datos nesecarios.

virus <- c("43740575","PP734539.1", "PP573039", "PP106453", "OR357774", "PP125282", "PP693357", "OL966962", "BS007901", "ON507015", "OZ027746", "PP754085", "OR529199", "ON115272", "OL869974", "PP038164", "OP012884", "OM287553", "PP107476", "OY778990", "ON239776")
virus_sequences<- ape::read.GenBank(virus)
write.dna(virus_sequences, file = "virus_sequences.fasta", format = "fasta")
fasta <- readDNAStringSet("virus_sequences.fasta", format = "fasta")

clean_sequences <- lapply(fasta, function(seq) {
  # Replace non-ATCG characters with nothing (remove them)
  clean_seq <- gsub("[^ATCGatcg]", "", as.character(seq))
  return(DNAString(clean_seq))
})
fasta <- DNAStringSet(clean_sequences, use.names = TRUE)



# Asignar nombres específicos a las subsecuencias
names(fasta) <- c("wuhan","usa", "china", "india", "francia", "alemania", "brasil", "corea_del_sur", "japon", "italia", "united_kingdom", "russia", "turkey", "spain", "australia", "vietnam", "argentina", "netherlands", "mexico", "dinamarca", "canada")

# Guardar las subsecuencias en un nuevo archivo FASTA
writeXStringSet(fasta, filepath = "virus_sequences.fasta", format = "fasta")
calculate_nucleotide_frequencies <- function(sequence) {
  bases <- as.character(sequence)
  freq <- table(strsplit(bases, ""))
  return(as.data.frame(freq))
}
fasta <- readDNAStringSet("virus_sequences.fasta", format = "fasta")
attributes(fasta)
calculate_gc_content <- function(sequence) {
  gc_count <- sum(letterFrequency(sequence, letters = c("G", "C"), as.prob = FALSE))
  total_bases <- length(as.character(sequence))
  return((gc_count / total_bases) * 100)
}
gc_contents <- sapply(fasta, calculate_gc_content)
gc_data <- data.frame(Virus = names(fasta), GC_Content = gc_contents)
# Define a custom color palette
custom_colors <- colorRampPalette(c("#FF0000", "#0000FF", "#00FF00", "#FF00FF", "#00FFFF", "#FFFF00", "#000000"))(length(unique(gc_data$Virus)))

# Use the custom palette in your plot
ggplot(gc_data, aes(x = Virus, y = GC_Content, fill = Virus)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = sprintf("%.2f%%", GC_Content)), vjust = -0.3, color = "black") +
  theme_minimal() +
  labs(title = "Porcentaje de Contenido de GC por Virus", x = "Virus", y = "Porcentaje de GC (%)") +
  scale_fill_manual(values = custom_colors)



# Aplicar la función a cada secuencia en el objeto fasta
nucleotide_frequencies <- lapply(fasta, calculate_nucleotide_frequencies)

# Crear un dataframe único para todas las secuencias
nucleotide_df <- do.call(rbind, nucleotide_frequencies)
nucleotide_df$Virus <- rep(names(fasta), sapply(nucleotide_frequencies, nrow))

# Ajustar nombres de columnas
names(nucleotide_df) <- c("Nucleotide", "Frequency", "Virus")

# Crear el gráfico de barras utilizando ggplot2
ggplot(nucleotide_df, aes(x = Nucleotide, y = Frequency, fill = Nucleotide)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Frequency), vjust = -0.3, color = "black") +  # Ajusta la posición vertical si es necesario
  facet_wrap(~ Virus, scales = "free_y") +
  theme_minimal() +
  labs(title = "Frecuencia de Nucleótidos por Virus",
       x = "Nucleótido", y = "Frecuencia") +
  scale_fill_brewer(palette = "Set1")

virus_seq_not_align <-
  Biostrings::readDNAStringSet("virus_sequences.fasta", format = "fasta")
str(virus_seq_not_align)
virus_seq_align_1 <- DECIPHER::AlignSeqs(virus_seq_not_align)

library(Biostrings)

# Assuming 'virus_protein' contains the aligned protein sequences as an AAStringSet
virus_protein <- Biostrings::translate(virus_seq_not_align)
names(virus_protein) <- names(virus_seq_not_align)

# Initialize an empty matrix to store distances
n <- length(virus_protein)
dist_matrix <- matrix(0, nrow = n, ncol = n)
rownames(dist_matrix) <- names(virus_protein)
colnames(dist_matrix) <- names(virus_protein)

# Calculate pairwise alignment scores
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    # The score is returned directly as a numeric value
    score <- pairwiseAlignment(virus_protein[[i]], virus_protein[[j]], substitutionMatrix = "BLOSUM62", gapOpening = -10, gapExtension = -1, scoreOnly = TRUE)
    # Convert score to a distance: here we use negative score for simplicity
    dist_matrix[i, j] <- -score
    dist_matrix[j, i] <- dist_matrix[i, j]  # Ensure the matrix is symmetric
  }
}

# Print the distance matrix
print(dist_matrix)

# Optional: Clustering or tree construction from the distance matrix
hc <- hclust(as.dist(dist_matrix), method = "complete")
plot(hc, main = "Hierarchical Clustering of Virus Proteins")

# Simple color palette using base R
colors <- colorRampPalette(c("white", "black"))(256)

# Apply this palette to your heatmap
heatmap_plot <- heatmap(as.matrix(dist_matrix), symm = TRUE, 
                        Rowv = as.dendrogram(hc), Colv = as.dendrogram(hc), 
                        col = colors, main = "Matriz de Distancia",
                        labRow = names(virus_protein), labCol = names(virus_protein),
                        cexRow = 0.5, cexCol = 0.5,  # Adjust text size for row and column labels
                        key = TRUE, keysize = 1.0)
print(heatmap_plot)

# Convert hierarchical clustering to dendrogram
dend <- as.dendrogram(hc)

# Plot the dendrogram vertically
plot(dend, main = "Hierarchical Clustering of Virus Proteins", 
     xlab = "Height", ylab = "Virus Proteins", 
     sub = "", ylim = c(0.5, length(virus_protein) + 0.5))
# Convert hclust to phylo
phylo_tree <- as.phylo(hc)
# Install necessary packages if not already installed
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
if (!requireNamespace("plotly", quietly = TRUE)) install.packages("plotly")  # For interactive functionality

# Load the packages
library(igraph)
library(ggraph)
library(plotly)

# Assuming 'aligned_sequences_dnabin' is available
graph <- graphMutations(aligned_sequences_dnabin, plot=FALSE)

# Assigning mutation types for illustration
E(graph)$mutation_type <- sample(c("substitution", "insertion", "deletion"), gsize(graph), replace=TRUE)

# Define a simplified color palette for the three mutation types
colors <- c("substitution" = "darkgreen", "insertion" = "darkorange", "deletion" = "black")

# Use ggraph to plot this igraph object with customized aesthetics
p <- ggraph(graph, layout = 'with_kk') +  # Kamada-Kawai layout for better spacing
  geom_edge_link(aes(edge_width = 0.1, color = mutation_type), alpha = 0.2, curved = FALSE) +  # Thinner, more transparent edges
  geom_node_point(size = 3, color = "red") +
  geom_node_text(aes(label = name), check_overlap = TRUE, color = "black", size = 3) +
  scale_edge_color_manual(values = colors) +
  theme_void() +
  labs(title = "Network Graph of Mutations", color = "Mutation Type")

# Print the static graph
print(p)

