import os

def run_tests(test_names):

    run_dir = "~/git/cea/build-dev/source"

    for test in test_names:
        # Execute the code on the input file
        print(f"Running {test}")
        subprocess.run(run_dir+"/cea"+f" {test}", shell=False, check=True)
        print()

    return

if __name__=="__main__":
    # test_names = ["example1", "example2", "example3",
    #               "example4","example5", "example6",
    #               "example7", "example8","example9",
    #               "example10","example11", "example12",
    #               "example13", "example14"]
    test_names = ["example1"]
    run_tests(test_names)