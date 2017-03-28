Here are some details on how we’re organizing the forum. In brief, the forum, will be a sequence of code demonstrations (with Q&A and discussion throughout), aimed to show novices how R can answer practical archaeological questions. The keywords here are ‘novice’ and ‘practical’, and the guidelines below are designed to ensure that we have maximum impact on archaeologists new to R (and that this impact is positive :)

I will be in touch with you all after the holidays to discuss the specific schedule of our 2-hour session.

Here are the guidelines for preparing for the forum “How to do archaeological science using R” (Session ID 2717)

- Prepare for 10 minutes of live-coding demonstration, and be willing to take questions throughout (but only take quick ones, defer complex questions until after your demo).     
- Prepare for others to code-along with you, using all kinds of different setups - aim for your code to work on vanilla R without custom shortcuts, config/profile modifications, etc. It’s ok to use CRAN packages.     
- Focus on the archaeological application, propose an archaeological question, and then use code to show an answer. We don’t want to get bogged down in programming esoterica     
- Many archaeologists already have the idea that R is hard, difficult to learn, and they can do everything they want in Excel. One of the goals of this forum is to challenge these notions. Image that our audience is somewhat skeptical, and that we need to persuade them that getting started with R is not hard, learning R can be fun, and that R is useful. This is not the time for exhaustive treatments of your topic. Instead, create a positive, upbeat and optimistic atmosphere for your demo, and make it a friendly introduction for our colleagues.     
- To achieve this goal, we must make some difficult trade-offs. We are making a deliberate decision here to hide some of the complexity and flexibility of R to instead highlight the accessibility and user-friendliness. Maybe in the future we can have a workshop on advanced R, but I think most archaeologists are a long way from this.     

## How to structure your demo:

- Start by explaining the archaeological question you will answer, give and brief context to motivate your example. It doesn’t have to be profound/novel/complicated, the only requirement here is that it is a realistic question and one you think (or know) that other archaeologists are interested in. It’s ok if you’ve already published the results.     
- Show data ingest from a local file or an online data source.     
- Show some basic data checking (e.g. str(), dplyr::glimpse(), View() )    
- Show some visualisation of the data, just a few simple plots to demonstrate good practices in Exploratory Data Analysis.    
- Show the packages that you use with library(xxx) or require(xxx)    
- We want to avoid magic tricks, ‘just do this’, or highly customised setups. To avoid these, will put long scripts into a function in a package (I can do this for you) so that we can use it as a one-liner in the forum. During the workshop, we’ll use the custom package and one-liner functions to abstract away some of the complexity. Think of the novices!    
- Omit extensive diagnostics. While these are important, they are better dealt with in an ‘intermediate’ or ‘advanced’ workshop. Here we’re aiming at complete novices, and we want to introduce them to common workflows with R.     
- Show the output in a form that can be easily adapted to other contexts (ie. pasted into a Word doc, sent by email). For example, make sure you output a plot as a PNG file, or output a table as a CSV file.     
- End by showing how the analysis has answered the archaeological question that you started with.     

## Code style requirements

- Write your code in a single .R or, ideally, .Rmd file (to knit to HTML, for speed and minimal dependencies). Prepare for big non-R packages such as TeX to be unavailable. Imaging that CRAN and github are all that we can use.     
- Use in-line comments in your code in plain English to help the audience follow along    
- Use highly readable code styling, especially blank lines between code chunks and spaces after commas and around operators (e.g. <-, =, +, etc.), indent liberally, and keep lines short to minimise side-to-side scrolling    

## Post-conference publication

Matt has been contacted by Springer to produce an edited volume based on the presentations in this forum. That has inspired us to also prepare an online book of your demonstrations as ‘Case studies in archaeological science using R’. This online book will also become the published book, but the online book will remain freely accessible and a lot easier to use. 

We will adapt the demos into short chapters of annotated tutorials (1000-2000 words) in an interactive book hosted on GitHub (see here for examples and more details: https://bookdown.org/). We’ll send more details on how to prepare your chapter, but it should be pretty simple. To prepare for this we will need to add some brief text to give background context and motivate and narrate our demos. Some citations will also be necessary. If you can’t make it to the SAA2017 forum, please consider making a submission to this edited volume. 

Don’t worry, this will be a lot less work than a regular book chapter. It is ok to reuse an analysis that has appeared elsewhere because the focus here is on a step-by-step walkthrough of the code.  
