# GCTA_doc
Document of GCTA

## Install NPM

## Install NPM packages
cd the cloned folder
```
npm install
```

## Revise the Markdown files
In md folder

We can use any markers from Markdown.

Some additional notes:

* Use ## as the header to show on the left navigation bar 

> Put {: .expand} to the end of ##, so the page will expand the contents not list the sub-titles when clicked;

* Use ### to show as sub-header title. 

> Put {: .notoc} to the end of ###, so this sub-title will not show in the navigation bar. 

We can not make these titles too long, because it will mess the navigation bar

* The hyperlinks are relative to build folder, we should put our resources,
such as gcta program into build folder after conver to HTML

We'd better preview it after editing. Some editor provides preview funtion. Pay attention to the contents, 
don't care about the content style, because it will be changed in our further steps further. 

## Convert to HTML
```
npm run build
```

## Prepare offline pdf

1. Goto build open the webpage, click the print in the browser, it will generate a nice PDF

2. Fix the PDF, as some blank pages, and no titiles

3. Put the PDF into build/static folder with the name gcta\_doc\_latest.pdf

4. Put the latest GCTA into build/static with the name gcta_latest.zip

## Create a zip of build folder, and send it to the web administrator.

> Note:
> Don't forget to put the new GCTA in the build folder
> Don't package the whole folder, only the **build** folder is needed for the deployment.


## Contact

If you have addtional problems, email zhilizheng@outlook.com for help.


