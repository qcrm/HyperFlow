CODE_SRC := ../../src/fvm_2d_cart/
INC_DIR := ../../include/

SOURCE := $(CODE_SRC)main.cpp \
		$(CODE_SRC)simulation/simulation/simulation.cpp \
		$(CODE_SRC)simulation/parser/parser.cpp \
		$(CODE_SRC)simulation/output/euler_data_output/euler_data_output.cpp \
		$(CODE_SRC)simulation/output/data_output/data_output.cpp \
		$(CODE_SRC)scheme/boundary/dirichlet/dirichlet_boundary.cpp \
		$(CODE_SRC)scheme/boundary/transmissive/transmissive_boundary.cpp \
		$(CODE_SRC)scheme/boundary/reflective/reflective_boundary.cpp \
		$(CODE_SRC)scheme/boundary/internal/internal_boundary.cpp \
		$(CODE_SRC)scheme/boundary/boundary/boundary.cpp \
		$(CODE_SRC)scheme/limiter/vanleer/vanleer.cpp \
		$(CODE_SRC)scheme/limiter/minbee/minbee.cpp \
		$(CODE_SRC)scheme/limiter/limiter/limiter.cpp \
		$(CODE_SRC)scheme/godunov/godunov.cpp \
		$(CODE_SRC)scheme/ode/forward_euler/forward_euler.cpp \
		$(CODE_SRC)scheme/ode/ode/ode.cpp \
		$(CODE_SRC)scheme/riemann/riemann_hllc/riemann_hllc.cpp \
		$(CODE_SRC)scheme/riemann/riemann/riemann.cpp \
		$(CODE_SRC)scheme/scheme/scheme.cpp \
		$(CODE_SRC)scheme/timestep/timestep.cpp \
		$(CODE_SRC)mesh/mesh/mesh.cpp \
		$(CODE_SRC)mesh/mesh_block/mesh_block.cpp \
		$(CODE_SRC)simulation/initcon/sphericalrp/sphericalrp.cpp \
		$(CODE_SRC)simulation/initcon/quadrp/quadrp.cpp \
		$(CODE_SRC)simulation/initcon/doublerp/doublerp.cpp \
		$(CODE_SRC)simulation/initcon/constant/constant.cpp \
		$(CODE_SRC)simulation/initcon/initcon/initcon.cpp \
		$(CODE_SRC)model/euler_ideal_gas/euler_ideal_gas.cpp \
		$(CODE_SRC)model/model/model.cpp \
		$(CODE_SRC)tensor/tensor.cpp

FLAGS := -O3 -std=c++11
INC := -I$(INC_DIR)

all: $(SOURCE)
	 g++ $(FLAGS) $(INC) -o ../../bin/hyperflow $(SOURCE)

clean:
	 rm -f $(CODE_SRC)/*.o
	 rm -f ../../bin/hyperflow
