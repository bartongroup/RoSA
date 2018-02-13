# This file is part of RoSA.
# 
# RoSA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# RoSA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RoSA  If not, see <http://www.gnu.org/licenses/>.

MAKE_ANNOT_PATH = file.path('..','python','rosa','make_annotation.py')
COUNT_SPLICED_PATH = file.path('..','python','rosa','antisense.py')


#' Build an antisense annotation
#' 
#' make_annotation creates an antisense annotation (as gtf) from 
#' a standard annotation (as gff or gtf), which can then be used to generate 
#' antisense read counts (input 2) via your favourite read counting tool 
#' (e.g. featureCounts). The annotation produced by make_annotation only 
#' contains antisense features and so cannot be used in place of a standard 
#' annotation.
#' *** Depends on python being installed! ***
#' @param annotation_file annotation file as gff or gtf
#' @param output_file output file name without file extension
#' @export
make_annotation <- function(annotation_file, output_file)
{
  system2('python', 
          args=c(MAKE_ANNOT_PATH,'-a',annotation_file,'-o',output_file),
          wait=TRUE)
}

#' Count spliced sense and antisense reads
#' 
#' count_spliced generates sense and antisense counts of reads at splice 
#' junctions. The script takes a standard annotation (as gtf/gff) 
#' and corresponding alignment (as bam) and outputs counts of spliced sense 
#' and antisense reads to a designated output file. An index file (.bai file) 
#' should also have been pre-generated and be in the same directory as the 
#' bam file. Because the script must process an entire bam file of reads, it 
#' is very slow. The script is set up to break the bam file into chunks and 
#' process each chunk separately using sambamba and some custom filtering. 
#' On a cluster with drmaa installed, the script will use drmaa to submit 
#' each chunk as a separate job. On a single machine, the script will spawn 
#' a new process to run each chunk separately. Once all of the jobs have run, 
#' the script collates the results to give a count of the spliced reads.
#' *** Depends on python and sambamba being installed! ***
#' @param annotation_file annotation file as gff or gtf
#' @param alignment_file alignment file as bam
#' @param output_file output file name without file extension
#' @export
count_spliced <- function(annotation_file, alignment_file, output_file)
{
  system2('python', 
          args=c(COUNT_SPLICED_PATH,'-a',annotation_file,'-l', alignment_file, '-o',output_file),
          wait=TRUE)
}