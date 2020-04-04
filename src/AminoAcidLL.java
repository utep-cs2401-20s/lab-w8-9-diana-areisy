class AminoAcidLL {
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL() {

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon) {

    aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    counts = new int[codons.length];
    incrementCodons(inCodon);
    next = null;

  }

  /********************************************************************************************/
  /* Helper method to keep the codon counter updated*/
  public void incrementCodons(String c) {
    for (int i = 0; i < this.codons.length; i++) {
      if (codons[i].equals(c)) {
        counts[i]++;
      }
    }
  }


  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops,
   * if not passes the task to the next node.
   * If there is no next node, add a new node to the list that would contain the codon.
   */
  private void addCodon(String inCodon) {
    if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
      incrementCodons(inCodon);
    } else {
      if (next != null) {
        next.addCodon(inCodon);
      } else {
        next = new AminoAcidLL(inCodon);
      }
    }
  }


  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount() {

    int sum = 0;

    for (int i = 0; i < counts.length; i++) {
      sum += counts[i];
    }
    return sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList) {
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
   *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList) {
    int diff = 0;
    for (int i = 0; i < codons.length; i++) {
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts.
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList) {
    if (!inList.isSorted()) {
      System.out.print("List is not sorted");
      return -1;
    }
    if (next == null && inList.next == null) {
      return totalDiff(inList);
    }
    if (inList.next == null) {
      return totalCount() + next.aminoAcidCompare(inList);
    }
    if (next == null) {
      return totalCount() + aminoAcidCompare(inList.next);
    }

    return totalDiff(inList) + (next.aminoAcidCompare(inList.next));
  }

  /********************************************************************************************/
  /* Same as above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList) {
    if (!inList.isSorted()) {
      System.out.print("List is not sorted");
      return -1;
    }
    if (next == null && inList.next == null) {
      return totalDiff(inList);
    }
    if (inList.next == null) {
      return totalCount() + next.codonCompare(inList);
    }
    if (next == null) {
      return totalCount() + codonCompare(inList.next);
    }

    return codonDiff(inList) + (next.codonCompare(inList.next));
  }

  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList() {
    //Base case, if the next is null
    if (next == null) {
      return new char[]{aminoAcid};
    }

    // If next is not null, used recursion to keep adding to the total
    char[] a = next.aminoAcidList();
    char[] total = new char[a.length + 1];
    total[0] = aminoAcid;
    for (int i = 0; i < a.length; i++) {
      total[i + 1] = a[i];
    }

    return total;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts() {
    // base, if next is null
    if (next == null) {
      return new int[]{totalCount()};
    }
    int[] a = next.aminoAcidCounts();
    int[] total = new int[a.length + 1];
    total[0] = totalCount();

    for (int i = 0; i < a.length; i++) {
      total[i + 1] = a[i];
    }
    return total;
  }

  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted() {
    //If it gets to the end of the list then it is sorted
    if (next == null) {
      return true;
    }
    //If the next node is greater than the current
    // node based on ascii value, then it is not sorted.
    else if (aminoAcid > next.aminoAcid) {
      return false;
    }
    //Used recursion to go through the linked list
    return next.isSorted();
  }

  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence) {
    AminoAcidLL list = new AminoAcidLL(inSequence.substring(0, 3));
    while (inSequence.length() > 3 && AminoAcidResources.getAminoAcidFromCodon(inSequence.substring(0, 3)) != '*') {
      inSequence = inSequence.substring(3);
      list.addCodon(inSequence.substring(0, 3));
    }
    return list;
  }

  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList) {
    //Base case if list has one or less
    if (inList == null || inList.next == null) {
      return inList;
    }

    //Used selection sort to sort the list
    char temp;
    for (AminoAcidLL i = inList; i.next != null; i = i.next) {
      for (AminoAcidLL j = i.next; j != null; j = j.next) {
        if (i.aminoAcid > j.aminoAcid) {
          temp = i.aminoAcid;
          i.aminoAcid = j.aminoAcid;
          j.aminoAcid = temp;
        }
      }
    }
      return inList;
  }
}
