version 1.0

workflow HelloWorld {
    input {
        File input_file
        Int scatter_num
    }

    scatter (i in range(scatter_num)) {
        call WriteGreeting as wg1 {
            input:
                input_file = input_file
        }
        call WriteGreeting as wg2 {
            input:
                input_file = wg1.output_greeting
        }
        call WriteGreeting as wg3 {
            input:
                input_file = wg2.output_greeting
        }
    }

    output {
        Array[File] output_file = wg3.output_greeting
    }
}

task WriteGreeting {
    input {
        File input_file
    }
    String ubuntu_docker = "ubuntu:20.04"
    command {
        echo "Hello World"
    }
    runtime {
        docker: ubuntu_docker
        disk: "50 GB"
        memory: "2000 MiB"
        cpu: 1
    }
    output {
        # Write output to standard out
        File output_greeting = stdout()
    }
}
