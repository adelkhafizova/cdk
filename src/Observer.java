import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.tools.SystemOutLoggingTool;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

public class Observer {
    Observer(String original, String trained) {
        original_data_file = original;
        trained_data_file = trained;
    }
    void initialize_algorithms(double p_value_threshold, int minimum_occurrence, double substructure_frequency) {
        algorithm_original = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency);
        algorithm_original.data_initialization(original_data_file, "Inhibitor", "Noninhibitor");
        algorithm_trained = new Algorithm(p_value_threshold, minimum_occurrence, substructure_frequency);
        algorithm_trained.data_initialization(trained_data_file, "Inhibitor", "Noninhibitor");
    }
    void analyze() {
        active_active = new HashMap<String, String>();
        active_inactive = new HashMap<String, String>();
        inactive_active = new HashMap<String, String>();
        inactive_inactive = new HashMap<String, String>();
        Iterator it = algorithm_trained.active.entrySet().iterator();
        float active = 0;
        float false_negative = 0;
        float false_positive = 0;
        float inactive = 0;
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                active_active.put((String) pair.getKey(), (String) pair.getValue());
                ++active;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                active_inactive.put((String)pair.getKey(), (String)pair.getValue());
                ++false_positive;
            }
        }
        it = algorithm_trained.inactive.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            if (algorithm_original.active.containsKey(pair.getKey())) {
                inactive_active.put((String)pair.getKey(), (String)pair.getValue());
                ++false_negative;
            }
            if (algorithm_original.inactive.containsKey(pair.getKey())) {
                inactive_inactive.put((String)pair.getKey(), (String)pair.getValue());
                ++inactive;
            }
        }
        //float active = (float)active_active.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_negative = (float)active_inactive.size()/(float)algorithm_original.final_signatures_active.size();
        //float false_positive = (float)inactive_active.size()/(float)algorithm_original.final_signatures_inactive.size();
        //float inactive = (float)inactive_inactive.size()/(float)algorithm_original.final_signatures_inactive.size();
        System.out.println(active);
        System.out.println(false_negative);
        System.out.println(false_positive);
        System.out.println(inactive);
    }
    void run() throws FileNotFoundException, UnsupportedEncodingException, CDKException{
        algorithm_original.run("");
        algorithm_trained.run("");
    }
    Algorithm algorithm_original;
    Algorithm algorithm_trained;
    String original_data_file;
    String trained_data_file;
    HashMap<String, String> active_active; //active originally and active in trained
    HashMap<String, String> active_inactive; //active originally but inactive in trained
    HashMap<String, String> inactive_inactive; // inactive originally and inactive in trained
    HashMap<String, String> inactive_active; //inactive originally and active in trained
}
