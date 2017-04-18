$(document).ready(function() {
    // if SVD is enabled (url), show SVD note
    if (RDC_calc) {
        console.log("SVD enabled!");
        $(".SVDnote").removeClass("hidden");
    }

    // If comp. measure is Q-VALUE, hide non-RDC inputs
    $(".mes").on( 'change', function() {
        console.log("MEASURE changed");
        var measures = document.querySelectorAll(".options");
        var measure;

        // get comp measure
        for (i = 0; i < measures.length; i++) {
            if (measures[i].checked) {
                measure = measures[i].value;
            }
        }

        if (measure == "q-value") {
            // show Q-VALUE note
            $(".QVALnote").removeClass("hidden");

            // hide non-RDC input fields
            $(".inputrange").each(function() {
                my_id    = $( this ).attr('id');
                my_range = $( this ).attr('id').split('_');
                my_type  = my_range[0];

                if (my_type != "RDC") {
                    // clear non-RDC inputs
                    $( this ).val(0);
                    $( this ).closest(".selection-table").addClass("hidden");
                }
            })
        // if compl. measure is other than Q-VALUE
        } else {
            // hide Q-VALUE note
            $(".QVALnote").addClass("hidden");

            // show non-RDC input fields
            $(".inputrange").each(function() {
                my_id    = $( this ).attr('id');
                my_range = $( this ).attr('id').split('_');
                my_type  = my_range[0];

                if (my_type != "RDC") {
                    $( this ).closest(".selection-table").removeClass("hidden");
                }
            })
        }
    });

    $("#bulk_RDC").click(function () {
        $(".inputrange").each(function() {
            my_id    = $( this ).attr('id');
            my_range = $( this ).attr('id').split('_');
            my_type  = my_range[0];

            if (my_type == "RDC") {
                $( this ).val("10");
                $('output[for="' + my_id + '"]').prop('value', "10");
            }
        })
        console.log("RDC bulk selection done.")
    });

    $("#bulk_S2").click(function () {
        $(".inputrange").each(function() {
            my_id    = $( this ).attr('id');
            my_range = $( this ).attr('id').split('_');
            my_type  = my_range[0];

            if (my_type == "S2") {
                $( this ).val("10");
                $('output[for="' + my_id + '"]').prop('value', "10");
            }
        })
        console.log("S2 bulk selection done.")
    });

    $("#bulk_Jcoup").click(function () {
        $(".inputrange").each(function() {
            my_id    = $( this ).attr('id');
            my_range = $( this ).attr('id').split('_');
            my_type  = my_range[0];

            if (my_type == "Jcoup") {
                $( this ).val("10");
                $('output[for="' + my_id + '"]').prop('value', "10");
            }
        })
        console.log("Jcoup bulk selection done.")
    });

    $("#bulk_CS").click(function () {
        $(".inputrange").each(function() {
            my_id    = $( this ).attr('id');
            my_range = $( this ).attr('id').split('_');
            my_type  = my_range[0];

            if (my_type == "CS") {
                $( this ).val("10");
                $('output[for="' + my_id + '"]').prop('value', "10");
            }
        })
        console.log("CS bulk selection done.")
    });


    $(".inputrange").on( 'change', function() {
        if (RDC_calc) {
            my_range = $( this ).attr('id').split('_');
            sel_type = my_range[0];
            sel_num  = my_range[1];
            sel_val  = $( this ).val();

            if (sel_type == "RDC") {
                $(".inputrange").each(function() {
                    my_id    = $( this ).attr('id');
                    my_range = $( this ).attr('id').split('_');
                    my_type  = my_range[0];
                    my_num   = my_range[1];

                    if (my_type == "RDC" && my_num == sel_num) {
                        $( this ).val(sel_val);
                        console.log('output[for="' + my_range + '"]');
                        $('output[for="' + my_id + '"]').prop('value', sel_val);
                    }
                })
            }
        }
    });

    $("#selection").click(function () {
        $(".involve-result").toggleClass("involve", 1000, "easeOutSine")
    });

    $('#startsel').click(function() {
        var ranges = document.querySelectorAll(".inputrange");
        var measures = document.querySelectorAll(".options");
        var min_size = $("#min_ens_size").val().trim();
        var max_size = $("#max_ens_size").val().trim();
        var overdrive = $("#overdrive").val().trim();
        var measure;
        var range_selected = false;
        var command = {};
        var i;

        // add selected measure to command
        for (i = 0; i < measures.length; i++) {
            if (measures[i].checked) {
                console.log("selected measure: " + measures[i].value);
                command["MEASURE"] = measures[i].value;
                measure = measures[i].value;
            }
        }

        if (min_size != "") {
            command["MIN_SIZE"] = min_size;
            console.log("minimum size is: " + min_size);
        }

        if (max_size != "") {
            command["MAX_SIZE"] = max_size;
            console.log("minimum size is: " + max_size);
        }

        if (overdrive != "") {
            command["OVERDRIVE"] = overdrive;
            console.log("overdrive is: " + overdrive);
        }

        // check if measure is selected
        if (measure == undefined) {
            alert("Please select compliance measure!");
            return false;
        }

        // add involved parameters to command
        for (i = 0; i < ranges.length; i++) {
            if (ranges[i].value != 0) {
                let myValue = String(ranges[i].id);
                command[myValue] = ranges[i].value;
                range_selected = true;
            }
        }

        // check if at least one parameter is involved in selection
        if (!range_selected) {
            alert("Please include at least one parameter!");
            return false;
        }

        // UPLAD HERE!
        console.log(command);

        let post_target = "/selection/" + $('#calculation_id').text();
        let selected_PDB_URL = "\"/media/" + $('#calculation_id').text() + "/selected.pdb\"";
        let my_id = $('#calculation_id').text();
        console.log("POST TARGET:", post_target, "ID", my_id);





        $.ajax({
            url: post_target,
            type: "POST",
            data: JSON.stringify(command),
            contentType: "application/json",
            complete: function(data, status){

                let sel_html = "\
                <h3 style=\"text-align: center;\">Selection results</h3>\
                <table class=\"files_table\">\
                <tr>\
                <td class=\"head-td\">" + data.responseJSON.measure + "</td>\
                <td>All models</td>\
                <td>Selected models: " + data.responseJSON.num_coordsets + "<a href=" + selected_PDB_URL + " download> (download)</a></td>\
                </tr>";


                for (let prop in data.responseJSON.values) {
                    console.log("obj." + prop + " = " + data.responseJSON.values[prop]);
                    sel_html += '<tr>';
                    sel_html += '<td class=\"head-td\" style=\"text-align: center;\">' + prop + '</td>';
                    sel_html += '<td style=\"text-align: center;\">' + data.responseJSON.values[prop].original + '</td>';
                    sel_html += '<td style=\"text-align: center;\">' + data.responseJSON.values[prop].selection + '</td>';
                    sel_html += '</tr>';
                }



                sel_html += "</table>";


                sel_html += "<div class=\"results\">\
                  <h4 class=\"table-source\">PCA projections of the <b>selected</b> ensemble</h4>";

                sel_html += '<table class=\"pca_table\"><tr>';
                sel_html += '<td class=\"pca_graph\"><a class=\"fancybox\" title=\"PCA mode 1-2\" href=\"/media/'
                         + my_id +'/pca_mode_12.svg\"><img width=\"320\" src=\"/media/'
                         + my_id +'/pca_mode_12.svg\"></td>';
                sel_html += '<td class=\"pca_graph\"><a class=\"fancybox\" title=\"PCA mode 2-3\" href=\"/media/'
                         + my_id +'/pca_mode_23.svg\"><img width=\"320\" src=\"/media/'
                         + my_id +'/pca_mode_23.svg\"></td>';
                sel_html += '<td class=\"pca_graph\"><a class=\"fancybox\" title=\"PCA mode 3-4\" href=\"/media/'
                         + my_id +'/pca_mode_34.svg\"><img width=\"320\" src=\"/media/'
                         + my_id +'/pca_mode_34.svg\"></td>';
                sel_html += '</tr></table>';

                sel_html += "</div>"


                sel_html += "<h3 style=\"text-align: center;\">Back-calculated data for the <b>original</b> ensemble</h3>";

                $(".selection-results").html(sel_html);
            }
        });


        return false;

    });
});
