.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('\n\nStego package\n\n')
        packageStartupMessage(
"                         .       .
                        / `.   .' \\
                .---.  <    > <    >  .---.
                |    \\  \\ - ~ ~ - /  /    |
                 ~-..-~             ~-..-~
             \\~~~\\.'                    `./~~~/
              \\__/                        \\__/
               /                  .-    .  \\
        _._ _.-    .-~ ~-.       /       }   \\/~~~/
    _.-'q  }~     /       }     {        ;    \\__/
   {'__,  /      (       /      {       /      `. ,~~|   .     .
    `''''='~~-.__(      /_      |      /- _      `..-'   \\\\   //
                / \\   =/  ~~--~~{    ./|    ~-.     `-..__\\\\_//_.-'
               {   \\  +\\         \\  =\\ (        ~ - . _ _ _..---~
               |  | {   }         \\   \\_\\
              '---.o___,'       .o___,'")
    }
}