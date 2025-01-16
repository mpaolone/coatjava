for dir in reconstruction/alert/src/main/java/org/jlab/rec/atof/*; do
    if [ -d "$dir/.git" ]; then
        git rm --cached "$dir"
        git config -f .gitmodules --remove-section submodule."$dir"
        rm -rf ".git/modules/$dir"
        git add "$dir"
    fi
done
git commit -m "Converted all submodules to regular directories"
git push origin development

