CC="/usr/bin/xcrun -sdk macosx10.13 clang"
                CXX="/usr/bin/xcrun -sdk macosx10.13 clang++"
                CFLAGS="-fno-common -arch x86_64 -mmacosx-version-min=10.9 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk -DMATLAB_MEX_FILE"
                CXXFLAGS="-fno-common -arch x86_64 -mmacosx-version-min=10.9 -fexceptions -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk -fobjc-arc -std=c++11 -stdlib=libc++ -DMATLAB_MEX_FILE"
                COPTIMFLAGS="-O2 -fwrapv -DNDEBUG"
                CXXOPTIMFLAGS="-O2 -fwrapv -DNDEBUG"
                CDEBUGFLAGS="-g"
                CXXDEBUGFLAGS="-g"
                LD="/usr/bin/xcrun -sdk macosx10.13 clang"
                LDXX="/usr/bin/xcrun -sdk macosx10.13 clang++"
                LDFLAGS="-Wl,-twolevel_namespace -undefined error -arch x86_64 -mmacosx-version-min=10.9 -Wl,-syslibroot,/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk -bundle  -Wl,-exported_symbols_list,"/Applications/MATLAB_R2018a.app/extern/lib/maci64/mexFunction.map" -L"/Applications/MATLAB_R2018a.app/bin/maci64" -lmx -lmex -lmat -lc++ -Wl,-exported_symbols_list,"/Applications/MATLAB_R2018a.app/extern/lib/maci64/mexFunction.map""
                LDDEBUGFLAGS="-g"
