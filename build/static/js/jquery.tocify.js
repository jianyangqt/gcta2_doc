!function(e){"use strict";!function(e,t,i,s){var o=".tocify-header",n=".tocify-subheader";e.widget("toc.tocify",{version:"1.9.0",options:{context:"body",tocID:"",ignoreSelector:null,showDIV:"",appendContents:"",threshLinks:1,expand:null,selectors:"h1, h2, h3",showAndHide:!0,showEffect:"slideDown",showEffectSpeed:"medium",hideEffect:"slideUp",hideEffectSpeed:"medium",smoothScroll:!0,smoothScrollSpeed:"medium",scrollTo:0,showAndHideOnScroll:!1,highlightOnScroll:!1,highlightOffset:40,theme:"bootstrap",extendPage:!0,extendPageOffset:100,history:!0,scrollHistory:!1,hashGenerator:"compact",highlightDefault:!0},_create:function(){var i=this;i.extendPageScroll=!0,i.items=[],i._generateToc(),i._addCSSClasses(),i.webkit=function(){for(var e in t)if(e&&-1!==e.toLowerCase().indexOf("webkit"))return!0;return!1}(),i._setEventHandlers(),e(t).on("hashchange load",function(){i._setActiveElement(!0),e("html, body").promise().done(function(){setTimeout(function(){i.extendPageScroll=!1},0)})})},_generateToc:function(){var t,i,s=this,o=s.options.ignoreSelector;if(t=-1!==this.options.selectors.indexOf(",")?e(this.options.context).find(this.options.selectors.replace(/ /g,"").substr(0,this.options.selectors.indexOf(","))):e(this.options.context).find(this.options.selectors.replace(/ /g,"")),!t.length)return void s.element.addClass("tocify-hide");s.element.addClass("tocify"),t.each(function(t){e(this).is(o)||(i=e("<ul/>",{id:"tocify-header"+t,class:"tocify-header"}).append(s._nestElements(e(this),t)),s.element.append(i),e(this).nextUntil(this.nodeName.toLowerCase()).each(function(){0===e(this).find(s.options.selectors).length?e(this).filter(s.options.selectors).each(function(){e(this).is(o)||s._appendSubheaders.call(this,s,i)}):e(this).find(s.options.selectors).each(function(){e(this).is(o)||s._appendSubheaders.call(this,s,i)})}))})},_setActiveElement:function(i){var s=this,o=t.location.hash.substring(1),n=s.element.find('li[data-unique="'+o+'"]');return o.length?(s.element.find("."+s.focusClass).removeClass(s.focusClass),n.addClass(s.focusClass),s.options.showAndHide&&n.click()):t.location.hash=e(s.options.tocID).children("ul").children("li").attr("data-unique"),s},_nestElements:function(t,i){var s,o,n;return s=e.grep(this.items,function(e){return e===t.text()}),s.length?this.items.push(t.text()+i):this.items.push(t.text()),n=this._generateHashValue(s,t,i),o=e("<li/>",{class:"tocify-item","data-unique":n}).append(e("<a/>",{text:t.text()})),t.before(e("<div/>",{name:n,"data-unique":n})),o},_generateHashValue:function(e,t,i){var s="",o=this.options.hashGenerator;if("pretty"===o){for(s=t.text().toLowerCase().replace(/\s/g,"-");s.indexOf("--")>-1;)s=s.replace(/--/g,"-");for(;s.indexOf(":-")>-1;)s=s.replace(/:-/g,"-")}else s="function"==typeof o?o(t.text(),t):t.text().replace(/\s/g,"");return e.length&&(s+=""+i),s},_appendSubheaders:function(t,i){var s=e(this).index(t.options.selectors),o=e(t.options.selectors).eq(s-1),a=+e(this).prop("tagName").charAt(1),l=+o.prop("tagName").charAt(1);a<l?t.element.find(n+"[data-tag="+a+"]").last().append(t._nestElements(e(this),s)):a===l?i.find(".tocify-item").last().after(t._nestElements(e(this),s)):i.find(".tocify-item").last().after(e("<ul/>",{class:"tocify-subheader","data-tag":a})).next(n).append(t._nestElements(e(this),s))},_setEventHandlers:function(){var s=this;this.element.on("click.tocify","li",function(i){if(s.options.history&&(t.location.hash=e(this).attr("data-unique")),s.element.find("."+s.focusClass).removeClass(s.focusClass),e(this).addClass(s.focusClass),s.options.showAndHide){var o=e('li[data-unique="'+e(this).attr("data-unique")+'"]');s._triggerShow(o)}s._scrollTo(e(this))}),this.element.find("li").on({"mouseenter.tocify":function(){e(this).addClass(s.hoverClass),e(this).css("cursor","pointer")},"mouseleave.tocify":function(){"bootstrap"!==s.options.theme&&e(this).removeClass(s.hoverClass)}}),(s.options.extendPage||s.options.highlightOnScroll||s.options.scrollHistory||s.options.showAndHideOnScroll)&&e(t).on("scroll.tocify",function(){e("html, body").promise().done(function(){var o,n,a,l,h=e(t).scrollTop(),r=e(t).height(),d=e(i).height(),c=e("body")[0].scrollHeight;if(s.options.extendPage&&(s.webkit&&h>=c-r-s.options.extendPageOffset||!s.webkit&&r+h>d-s.options.extendPageOffset)&&!e(".tocify-extend-page").length){if(n=e('div[data-unique="'+e(".tocify-item").last().attr("data-unique")+'"]'),!n.length)return;a=n.offset().top,e(s.options.context).append(e("<div />",{class:"tocify-extend-page",height:Math.abs(a-h)+"px","data-unique":"tocify-extend-page"})),s.extendPageScroll&&(l=s.element.find("li.active"),s._scrollTo(e('div[data-unique="'+l.attr("data-unique")+'"]')))}setTimeout(function(){var i,n=null,a=null,l=e(s.options.context).find("div[data-unique]");l.each(function(t){var i=Math.abs((e(this).next().length?e(this).next():e(this)).offset().top-h-s.options.highlightOffset);if(!(null==n||i<n))return!1;n=i,a=t}),i=e(l[a]).attr("data-unique"),o=e('li[data-unique="'+i+'"]'),s.options.highlightOnScroll&&o.length&&(s.element.find("."+s.focusClass).removeClass(s.focusClass),o.addClass(s.focusClass)),s.options.scrollHistory&&t.location.hash!=="#"+i&&t.location.replace("#"+i),s.options.showAndHideOnScroll&&s.options.showAndHide&&s._triggerShow(o,!0)},0)})})},show:function(t,i){var s=this;if(!t.is(":visible"))switch(t.find(n).length||t.parent().is(o)||t.parent().is(":visible")?t.children(n).length||t.parent().is(o)||(t=t.closest(n)):t=t.parents(n).add(t),s.options.showEffect){case"none":t.show();break;case"show":t.show(s.options.showEffectSpeed);break;case"slideDown":t.slideDown(s.options.showEffectSpeed);break;case"fadeIn":t.fadeIn(s.options.showEffectSpeed);break;default:t.show()}return t.parent().is(o)?s.hide(e(n).not(t)):s.hide(e(n).not(t.closest(o).find(n).not(t.siblings()))),s},hide:function(e){var t=this;switch(t.options.hideEffect){case"none":e.hide();break;case"hide":e.hide(t.options.hideEffectSpeed);break;case"slideUp":e.slideUp(t.options.hideEffectSpeed);break;case"fadeOut":e.fadeOut(t.options.hideEffectSpeed);break;default:e.hide()}return t},_triggerShow:function(e,t){var i=this;return e.parent().is(o)||e.next().is(n)?i.show(e.next(n),t):e.parent().is(n)&&i.show(e.parent(),t),i},_addCSSClasses:function(){return"jqueryui"===this.options.theme?(this.focusClass="ui-state-default",this.hoverClass="ui-state-hover",this.element.addClass("ui-widget").find(".toc-title").addClass("ui-widget-header").end().find("li").addClass("ui-widget-content")):"bootstrap"===this.options.theme?(this.element.find(o+","+n).addClass("nav nav-list"),this.focusClass="active"):(this.focusClass="tocify-focus",this.hoverClass="tocify-hover"),this},setOption:function(){e.Widget.prototype._setOption.apply(this,arguments)},setOptions:function(){e.Widget.prototype._setOptions.apply(this,arguments)},_scrollTo:function(t){var i=this,s=i.options.smoothScroll||0,o=(i.options.scrollTo,e('div[data-unique="'+t.attr("data-unique")+'"]')),n="";if(-1!==t.parent().attr("class").indexOf("tocify-header")){var a=t.parent(),l=a.children("ul").children("li");if(l.size()>i.options.threshLinks&&!o.next().is(i.options.expand))n+="<h2>"+a.children("li").children("a").text()+"</h2><ul>",l.each(function(t){n+="<li><a href='#"+e(this).attr("data-unique")+"' class='link_list'>"+e(this).text()+"</a></li>"}),n+="</ul>";else{var h,r=a.next().children("li");h=0!=r.size()?o.nextUntil('div[data-unique="'+r.attr("data-unique")+'"]').clone():o.nextAll().clone(),e(i.options.showDIV).html(h),n+=e(i.options.showDIV).html()}}else for(var d=o;;){var r=d.next();if(r.attr("data-unique"))break;if(0==r.size())break;n+=r[0].outerHTML,d=r}return n+=i.options.appendContents,e(i.options.showDIV).html(n),o.length?(e("html, body").promise().done(function(){e("html, body").animate({scrollTop:"0px"},{duration:s})}),i):i}})}(window.jQuery,window,document)}();