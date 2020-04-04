import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;

public class AminoAcidLLTester {

        /*
         * Testing aminoAcidList() method:
         * Making sure methods works properly by
         * recognizing the amino acids from a string
         * PASSED
         */
        @Test
        public void aminoAcidList() {
            String sequence = "AGGGAGGACUCA";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            char[] exp = {'R','E','D','S'};
            char[] result = list.aminoAcidList();
            assertArrayEquals(exp,result);
        }
        /*
         * Testing aminoAcidList() method:
         * PASSED
         */
        @Test
        public void aminoAcidList2() {
            String sequence = "GUGAUCCUGGAA";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            char[] exp = {'V','I','L','E'};
            char[] result = list.aminoAcidList();
            assertArrayEquals(exp,result);
        }
        /*
         * Testing aminoAcidCounts() method:
         * PASSED
         */
        @Test
        public void aminoAcidCounts1() {
            int[] expected = {2, 2, 1};
            String testSequence = "GCGGCGCAUCAUUGU";
            AminoAcidLL test = AminoAcidLL.createFromRNASequence(testSequence);
            assertArrayEquals(expected, test.aminoAcidCounts());
        }
        /*
         * Testing aminoAcidCounts() method:
         * PASSED
         * Set it up a different way and it still worked
         */
        @Test
        public void aminoAcidCounts2() {
            String sequence = "AAAAAAAGCCAACAACAAACA";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            int[] result = list.aminoAcidCounts();
            int[] exp = {2,1,3,1};
            assertArrayEquals(exp,result);
        }
        /*
         * Testing aminoAcidCompare() method:
         * PASSED
         * Since the amino acids are the same, the difference should be 0
         */
        @Test
        public void aminoAcidCompare() {
            AminoAcidLL list = AminoAcidLL.createFromRNASequence("CAUUUGAUU");
            AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("CAUUUGAUU");
            list = AminoAcidLL.sort(list);
            list2 = AminoAcidLL.sort(list2);
            assertEquals(0, list.codonCompare(list2));
        }
        /*
         * Testing codonCompare() method:
         * PASSED
         * Difference should be one since the second list has two codons and list has only one.
         */
        @Test
        public void codonCompare() {
            AminoAcidLL list = AminoAcidLL.createFromRNASequence("UUU");
            AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("UUUGCC");
            list = AminoAcidLL.sort(list);
            list2 = AminoAcidLL.sort(list2);
            assertEquals(1, list.codonCompare(list2));
        }
        /*
         * Testing isSorted() method:
         * PASSED
         * If the method returns true, then
         * it went through the whole string which means it is sorted.
         */
        @Test
        public void IsSorted1() {
            String sequence = "GCUUGUUUUUAC";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            assertEquals(true,list.isSorted());
        }
        /*
         * Testing isSorted() method:
         *PASSED
         */
        @Test
        public void IsSorted2() {
            String sequence = "GGCGCCACCCAA";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            assertEquals(false,list.isSorted());
        }
        /*
         * Testing sort() method:
         * Method should take the string, find the amino acids and place them in order.
         * PASSED
         */
        @Test
        public void Sort1() {
            AminoAcidLL string = AminoAcidLL.createFromRNASequence("ACAGCAGUUGAC"); //TAVD
            string = AminoAcidLL.sort(string);
            //Sorted:ADTV
            assertEquals(true, string.isSorted());

        }
        /*
         * Testing sort() method:
         * PASSED
         */
        @Test
        public void Sort2() {
            AminoAcidLL test = AminoAcidLL.createFromRNASequence("UUUUGACAU"); //FSTOPH
            test = AminoAcidLL.sort(test);
            //FHSTOP
            assertEquals(true, test.isSorted());
        }

        /*
         * Testing createFromRNASequence() method:
         * PASSED
         */
        @Test
        public void createFromRNASequence1() {
            String sequence = "AAAGAGAAUACU";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            AminoAcidLL temp = list;
            while (temp != null) {
                System.out.print(temp.aminoAcid + "   ");
                temp = temp.next;
                //It should print K E N T
            }
        }
        /*
         * Testing createFromRNASequence() method:
         * PASSED
         */
        @Test
        public void createFromRNASequence2() {
            String sequence = "CCCGCCAGAAAGUCA";
            AminoAcidLL list = AminoAcidLL.createFromRNASequence(sequence);
            AminoAcidLL temp = list;
            while (temp != null) {
                System.out.print(temp.aminoAcid + "   ");
                temp = temp.next;
                //It should print P A R K S
            }
        }
    }
