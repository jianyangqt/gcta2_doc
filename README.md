# GCTA_doc
Document of GCTA

## Install node.js and npm 

[Click here for instructions](https://nodejs.org/en/)

For Windows users, additional packages are required, click this link for [GNU tool chains under Windows](http://www.mingw.org/wiki/msys) 

## Install the embeded npm packages
Goto the cloned folder in terminal and run:
```
npm install
```

## Revise the Markdown files
Revise the Markdown file in the **md** folder. We can use any markers from Markdown.

Some additional notes:

* Use ## before the header to show on the left navigation bar 

> Put {: .expand} to the end of ##-prefixed header, if we want to expand the whole contents in that header title other than a list of 
the sub-titles when clicked in the navigation bar (e.g. 1overview.md);

* Use ### to show as sub-header title. 

> Put {: .notoc} to the end of ###, if we don't like this sub-title to show in the navigation bar. 

We can not make these titles too long, because it will mess the navigation bar

* The hyperlinks are relative to build folder, we should put our resources,
such as gcta program into build folder after conver to HTML

We'd better preview it after editing. Some editor provides preview funtion. Pay attention to the contents, 
don't care about the style, because it will be changed in our further steps. 

*Note: the Markdown files will be joined into a single file. Pay attention to the order of file name, thus it have 1 2 3 prefix.* 

## Convert to HTML
```
npm run build
```

## Prepare offline pdf

1. Goto **build** folder, open the webpage (index.html), click print in the web browser, it will generate a nice PDF

2. Fix the PDF, as there may have some blank pages, and no titiles

3. Put the PDF into build/static folder with the name gcta\_doc\_latest.pdf

4. Put the latest GCTA into build/static with the name gcta_latest.zip

*Note: this offline documents are for some users that have very old browsers, it is not prefered thus invisible to ordinary users*

## Create a zip of build folder

Send the zip file to the web administrator.

> Note:  
> Don't forget to put the new GCTA in the build folder  
> Don't package the whole folder, only the **build** folder is needed for the deployment.


## Help

If you have addtional problems, email zhilizheng@uq.edu.au (or zhilizheng@outlook.com) for help.

If we want to change the style of the web site, we need change the **template**, don't revise the builded web page directly.
