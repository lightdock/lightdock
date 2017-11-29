# How to contribute

Contributing to this project is more than welcome. The easiest way is to write your custom scoring functions and pull them for request.

## Custom scoring functions

You can write your own scoring function starting from the provided [driver.py](lightdock/scoring/template/driver.py) template.
The easiest way is to copy the template folder:

```
cp -R template/ mynewscoring
```

In order to call this new scoring function, just use the ```-s``` flag:

```
lightdock 2UUY_rec.pdb 2UUY_lig.pdb 1 10 5 -s mynewscoring
```
