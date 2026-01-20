
use chrono::Local;
use std::io::Write;

pub fn print_ascii_art( stdout: & mut Box<dyn Write>, low_memory: bool) {


    let start_time = Local::now();



    let ascii_art: &str = r#"
                                     ,---,.                            ,---,     
       ,---.                       ,'  .'  \  ,--,                   ,--.' |     
      /__./|   ,---.             ,---.' .' |,--.'|    __  ,-.        |  |  :     
 ,---.;  ; |  '   ,'\ ,--,  ,--, |   |  |: ||  |,   ,' ,'/ /|        :  :  :     
/___/ \  | | /   /   ||'. \/ .`| :   :  :  /`--'_   '  | |' | ,---.  :  |  |,--. 
\   ;  \ ' |.   ; ,. :'  \/  / ; :   |    ; ,' ,'|  |  |   ,'/     \ |  :  '   | 
 \   \  \: |'   | |: : \  \.' /  |   :     \'  | |  '  :  / /    / ' |  |   /' : 
  ;   \  ' .'   | .; :  \  ;  ;  |   |   . ||  | :  |  | ' .    ' /  '  :  | | | 
   \   \   '|   :    | / \  \  \ '   :  '; |'  : |__;  : | '   ; :__ |  |  ' | : 
    \   `  ; \   \  /./__;   ;  \|   |  | ; |  | '.'|  , ; '   | '.'||  :  :_:,' 
     :   \ |  `----' |   :/\  \ ;|   :   /  ;  :    ;---'  |   :    :|  | ,'     
      '---"          `---'  `--` |   | ,'   |  ,   /        \   \  / `--''       
                                 `----'      ---`-'          `----'              

                                 "#;
    
    let mut low_mem_msg = "";

    if low_memory {
        low_mem_msg = "Low memory mode\n"
    }

    writeln!(
        stdout,
        "{}\nCode Written by: Steven Pak\n{}Start date-time: {}", 
        ascii_art, 
        low_mem_msg,
        start_time.format("%Y-%m-%d %H:%M:%S")
    ).unwrap();
}