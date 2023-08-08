version 1.0

workflow HelloWorld {
    input {
        File input_file
    }

    call WriteGreeting {
        input:
            input_file = input_file
    }

    output {
        File output_file = WriteGreeting.output_greeting
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
        disks: "local-disk ${disk_size_gb} HDD"
        disk: "50 GB"
        memory: "2000 MiB"
        cpu: 1
    }
    output {
        # Write output to standard out
        File output_greeting = stdout()
    }
}
