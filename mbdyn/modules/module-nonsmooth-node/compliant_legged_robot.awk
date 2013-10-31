BEGIN {
	node_add("base_1", 0, -0.25, -0.15, 0, "hide");
	node_add("base_2", 0, 0.25, -0.15, 0, "hide");
	node_add("base_3", 0, 0.25, 0.15, 0, "hide");
	node_add("base_4", 0, -0.25, 0.15, 0, "hide");

	sideprop_add("base", 2);
	nodes[1] = "base_1";
	nodes[2] = "base_2";
	nodes[3] = "base_3";
	nodes[4] = "base_4";
	side_add("base", 4, nodes, "base");

	node_add("body_l_1", 100, -.075, -.04, -.03, "hide");
	node_add("body_l_2", 100, .075, -.04, -.03, "hide");
	node_add("body_l_3", 100, .075, .04, -.03, "hide");
	node_add("body_l_4", 100, -.075, .04, -.03, "hide");

	node_add("body_u_1", 100, -.075, -.04, .03, "hide");
	node_add("body_u_2", 100, .075, -.04, .03, "hide");
	node_add("body_u_3", 100, .075, .04, .03, "hide");
	node_add("body_u_4", 100, -.075, .04, .03, "hide");

	sideprop_add("body", 4);

	nodes[1] = "body_l_4";
	nodes[2] = "body_l_3";
	nodes[3] = "body_l_2";
	nodes[4] = "body_l_1";
	side_add("body_l", 4, nodes, "body");

	nodes[1] = "body_l_1";
	nodes[2] = "body_l_2";
	nodes[3] = "body_u_2";
	nodes[4] = "body_u_1";
	side_add("body_s", 4, nodes, "body");

	nodes[1] = "body_l_2";
	nodes[2] = "body_l_3";
	nodes[3] = "body_u_3";
	nodes[4] = "body_u_2";
	side_add("body_e", 4, nodes, "body");

	nodes[1] = "body_l_3";
	nodes[2] = "body_l_4";
	nodes[3] = "body_u_4";
	nodes[4] = "body_u_3";
	side_add("body_n", 4, nodes, "body");

	nodes[1] = "body_l_4";
	nodes[2] = "body_l_1";
	nodes[3] = "body_u_1";
	nodes[4] = "body_u_4";
	side_add("body_w", 4, nodes, "body");

	nodes[1] = "body_u_1";
	nodes[2] = "body_u_2";
	nodes[3] = "body_u_3";
	nodes[4] = "body_u_4";
	side_add("body_u", 4, nodes, "body");
}
