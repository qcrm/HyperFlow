CODE_SRC := ../../src/fvm_2d_uns/
INC_DIR := ../../include/

SOURCE := $(CODE_SRC)main.cpp \
		$(CODE_SRC)simulation/simulation/simulation.cpp \
		$(CODE_SRC)simulation/output/output.cpp \
		$(CODE_SRC)scheme/godunov/godunov.cpp \
		$(CODE_SRC)scheme/ode/forward_euler/forward_euler.cpp \
		$(CODE_SRC)scheme/ode/ode/ode.cpp \
		$(CODE_SRC)scheme/riemann/riemann_hllc/riemann_hllc.cpp \
		$(CODE_SRC)scheme/riemann/riemann/riemann.cpp \
		$(CODE_SRC)scheme/scheme/scheme.cpp \
		$(CODE_SRC)scheme/timestep/timestep.cpp \
		$(CODE_SRC)mesh/mesh_structured/mesh_structured.cpp \
		$(CODE_SRC)mesh/mesh/mesh.cpp \
		$(CODE_SRC)simulation/initcon/sphericalrp/sphericalrp.cpp \
		$(CODE_SRC)simulation/initcon/constant/constant.cpp \
		$(CODE_SRC)simulation/initcon/initcon/initcon.cpp \
		$(CODE_SRC)model/euler_ideal_gas/euler_ideal_gas.cpp \
		$(CODE_SRC)model/model/model.cpp \
		$(CODE_SRC)mesh/cell/cell.cpp \
		$(CODE_SRC)geometry/polygon/polygon.cpp \
		$(CODE_SRC)geometry/edge/edge.cpp \
		$(CODE_SRC)../share/tensor/tensor.cpp

FLAGS := -O3 -std=c++11
INC := -I$(INC_DIR)

all: $(SOURCE)
	 g++ $(FLAGS) $(INC) -o ../../bin/hyperflow2duns $(SOURCE)

clean:
	 rm -f $(CODE_SRC)/*.o
	 rm -f ../../bin/hyperflow2duns
