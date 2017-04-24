var jsdom = require("jsdom");
var fs = require("fs");
var copydir = require("copy-dir");

fs.readFile("whole.html", "utf8", function(error, data) {
    jsdom.env(data, [], function(errors, window) {
        var $ = require('jquery')(window);
        $("h3").before("<hr/>");
        $("p").each(function() {
            var content = $(this).html();
            var re = /--[a-zA-Z0-9\-]+/g;
            var re_content = content.replace(re, function(str) { return '<span class="kwd">' + str + '</span>'; });
            $(this).html(re_content);
        });


        $("h3:contains('{: '), h2:contains('{: ')").each(function() {
            var h_text = $(this).text();
            console.log("F :" + h_text + ".");
            var re_class = /\s*(.*?)\{:\s*\.(.*?)\}\s*/g;
            var match = re_class.exec(h_text);
            console.log(match);
            $(this).text(match[1]);
            $(this).addClass(match[2]);
        });

        fs.writeFile("whole_moded.html", window.document.documentElement.outerHTML, function(error) { if (error) throw error; });

        fs.readFile("./template/index.html", "utf8", function(error2, data2) {
            jsdom.env(data2, [], function(errors2, window2) {
                var $2 = require("jquery")(window2);
                $2('#content').html($('body'));

                fs.writeFile("./build/index.html", window2.document.documentElement.outerHTML, function(error) { if (error) throw error; });
                copydir.sync("./template/static", "./build/static");
                console.log("Done! See build folder");
            });
        });
    });
});