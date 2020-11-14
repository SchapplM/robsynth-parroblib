% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:20
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3PRRRR8V2G4P2A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(8,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4P2A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4P2A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4P2A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4P2A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4P2A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:20:21
% EndTime: 2020-08-06 18:20:21
% DurationCPUTime: 0.48s
% Computational Cost: add. (378->109), mult. (804->235), div. (9->3), fcn. (885->34), ass. (0->115)
t71 = legFrame(3,3);
t46 = sin(t71);
t52 = cos(t71);
t67 = sin(pkin(8));
t69 = cos(pkin(8));
t28 = t46 * t69 + t67 * t52;
t68 = sin(pkin(4));
t117 = t28 * t68;
t72 = legFrame(2,3);
t47 = sin(t72);
t53 = cos(t72);
t29 = t47 * t69 + t67 * t53;
t116 = t29 * t68;
t73 = legFrame(1,3);
t48 = sin(t73);
t54 = cos(t73);
t30 = t48 * t69 + t67 * t54;
t115 = t30 * t68;
t86 = cos(qJ(3,3));
t114 = pkin(3) * t86 ^ 2;
t88 = cos(qJ(3,2));
t113 = pkin(3) * t88 ^ 2;
t90 = cos(qJ(3,1));
t112 = pkin(3) * t90 ^ 2;
t111 = pkin(3) * t68;
t80 = sin(qJ(3,3));
t110 = t80 * pkin(2);
t82 = sin(qJ(3,2));
t109 = t82 * pkin(2);
t84 = sin(qJ(3,1));
t108 = t84 * pkin(2);
t31 = -t67 * t46 + t52 * t69;
t107 = t31 * t68;
t32 = -t67 * t47 + t53 * t69;
t106 = t32 * t68;
t33 = -t67 * t48 + t54 * t69;
t105 = t33 * t68;
t70 = cos(pkin(4));
t104 = t70 * t80;
t81 = sin(qJ(2,3));
t103 = t70 * t81;
t102 = t70 * t82;
t83 = sin(qJ(2,2));
t101 = t70 * t83;
t100 = t70 * t84;
t85 = sin(qJ(2,1));
t99 = t70 * t85;
t98 = t81 * t68;
t97 = t83 * t68;
t96 = t85 * t68;
t87 = cos(qJ(2,3));
t92 = pkin(7) + pkin(6);
t40 = pkin(2) * t81 - t92 * t87;
t95 = t80 * t111 - t40 * t70;
t89 = cos(qJ(2,2));
t41 = pkin(2) * t83 - t92 * t89;
t94 = t82 * t111 - t41 * t70;
t91 = cos(qJ(2,1));
t42 = pkin(2) * t85 - t92 * t91;
t93 = t84 * t111 - t42 * t70;
t79 = legFrame(1,2);
t78 = legFrame(2,2);
t77 = legFrame(3,2);
t76 = legFrame(1,1);
t75 = legFrame(2,1);
t74 = legFrame(3,1);
t63 = cos(t79);
t62 = cos(t78);
t61 = cos(t77);
t60 = sin(t79);
t59 = sin(t78);
t58 = sin(t77);
t57 = cos(t76);
t56 = cos(t75);
t55 = cos(t74);
t51 = sin(t76);
t50 = sin(t75);
t49 = sin(t74);
t45 = pkin(2) * t91 + t85 * t92;
t44 = pkin(2) * t89 + t83 * t92;
t43 = pkin(2) * t87 + t81 * t92;
t39 = t67 * t91 + t69 * t99;
t38 = t69 * t101 + t67 * t89;
t37 = t69 * t103 + t67 * t87;
t36 = t67 * t99 - t69 * t91;
t35 = t67 * t101 - t69 * t89;
t34 = t67 * t103 - t69 * t87;
t27 = pkin(3) * t100 + t68 * t42;
t26 = pkin(3) * t102 + t68 * t41;
t25 = pkin(3) * t104 + t68 * t40;
t24 = -t60 * t115 + t70 * t63;
t23 = -t59 * t116 + t70 * t62;
t22 = -t58 * t117 + t70 * t61;
t21 = t67 * t45 - t93 * t69;
t20 = t67 * t44 - t94 * t69;
t19 = t67 * t43 - t95 * t69;
t18 = -t45 * t69 - t93 * t67;
t17 = -t44 * t69 - t94 * t67;
t16 = -t43 * t69 - t95 * t67;
t15 = -t48 * t36 + t39 * t54;
t14 = -t47 * t35 + t38 * t53;
t13 = -t46 * t34 + t37 * t52;
t12 = 0.1e1 / (pkin(2) * t100 + t96 * t112 + t27 * t90);
t11 = 0.1e1 / (pkin(2) * t102 + t97 * t113 + t26 * t88);
t10 = 0.1e1 / (pkin(2) * t104 + t98 * t114 + t25 * t86);
t9 = t63 * t96 + (t36 * t54 + t39 * t48) * t60;
t8 = t62 * t97 + (t35 * t53 + t38 * t47) * t59;
t7 = t61 * t98 + (t34 * t52 + t37 * t46) * t58;
t6 = -t18 * t48 + t21 * t54;
t5 = -t17 * t47 + t20 * t53;
t4 = -t16 * t46 + t19 * t52;
t3 = t27 * t63 + (t18 * t54 + t21 * t48) * t60;
t2 = t26 * t62 + (t17 * t53 + t20 * t47) * t59;
t1 = t25 * t61 + (t16 * t52 + t19 * t46) * t58;
t64 = [(-((t30 * t99 - t91 * t33) * t63 - t60 * t96) * t112 + ((t93 * t30 + t33 * t45) * t63 + t27 * t60) * t90 + (t63 * t115 + t70 * t60) * t108) * t12, (-(-t57 * t15 + t9 * t51) * t112 + (-t3 * t51 + t6 * t57) * t90 - (t57 * t105 + t24 * t51) * t108) * t12, ((t15 * t51 + t9 * t57) * t112 + (t3 * t57 + t6 * t51) * t90 + (-t51 * t105 + t24 * t57) * t108) * t12; (-((t29 * t101 - t89 * t32) * t62 - t59 * t97) * t113 + ((t94 * t29 + t32 * t44) * t62 + t26 * t59) * t88 + (t62 * t116 + t70 * t59) * t109) * t11, (-(-t56 * t14 + t8 * t50) * t113 + (-t2 * t50 + t5 * t56) * t88 - (t56 * t106 + t23 * t50) * t109) * t11, ((t14 * t50 + t8 * t56) * t113 + (t2 * t56 + t5 * t50) * t88 + (-t50 * t106 + t23 * t56) * t109) * t11; (-((t28 * t103 - t87 * t31) * t61 - t58 * t98) * t114 + ((t95 * t28 + t31 * t43) * t61 + t25 * t58) * t86 + (t61 * t117 + t70 * t58) * t110) * t10, (-(-t55 * t13 + t7 * t49) * t114 + (-t1 * t49 + t4 * t55) * t86 - (t55 * t107 + t22 * t49) * t110) * t10, ((t13 * t49 + t7 * t55) * t114 + (t1 * t55 + t4 * t49) * t86 + (-t49 * t107 + t22 * t55) * t110) * t10;];
Jinv  = t64;
