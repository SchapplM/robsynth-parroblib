% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [4x4]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:33
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P4PRRRR8V2G3P1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(8,1),zeros(4,3),zeros(4,3)}
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G3P1A1_Jinv: qJ has to be [3x4] (double)');
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G3P1A1_Jinv: xP has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G3P1A1_Jinv: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G3P1A1_Jinv: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G3P1A1_Jinv: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:33:04
% EndTime: 2020-08-07 11:33:05
% DurationCPUTime: 0.92s
% Computational Cost: add. (388->157), mult. (864->331), div. (20->8), fcn. (854->30), ass. (0->134)
t59 = sin(qJ(2,4));
t61 = cos(qJ(2,4));
t78 = pkin(7) + pkin(6);
t155 = pkin(2) * t59 - t61 * t78;
t67 = sin(qJ(2,3));
t73 = cos(qJ(2,3));
t154 = pkin(2) * t67 - t73 * t78;
t69 = sin(qJ(2,2));
t75 = cos(qJ(2,2));
t153 = pkin(2) * t69 - t75 * t78;
t74 = cos(qJ(3,2));
t143 = pkin(3) * t74 ^ 2;
t110 = t69 * t143;
t152 = t153 * t74 + t110;
t72 = cos(qJ(3,3));
t144 = pkin(3) * t72 ^ 2;
t111 = t67 * t144;
t151 = t154 * t72 + t111;
t60 = cos(qJ(3,4));
t145 = pkin(3) * t60 ^ 2;
t112 = t59 * t145;
t150 = t155 * t60 + t112;
t58 = sin(qJ(3,4));
t149 = pkin(2) * t58;
t66 = sin(qJ(3,3));
t148 = pkin(2) * t66;
t68 = sin(qJ(3,2));
t147 = pkin(2) * t68;
t70 = sin(qJ(3,1));
t146 = pkin(2) * t70;
t76 = cos(qJ(3,1));
t142 = pkin(3) * t76 ^ 2;
t55 = sin(pkin(4));
t141 = pkin(3) * t55;
t56 = cos(pkin(8));
t140 = t56 * pkin(2);
t139 = t60 * pkin(3);
t138 = t72 * pkin(3);
t137 = t74 * pkin(3);
t54 = sin(pkin(8));
t136 = t54 * t55;
t57 = cos(pkin(4));
t135 = t54 * t57;
t134 = t55 * t59;
t133 = t55 * t67;
t132 = t55 * t69;
t131 = t55 * t70;
t71 = sin(qJ(2,1));
t130 = t55 * t71;
t129 = t56 * t78;
t128 = t57 * t59;
t127 = t57 * t67;
t126 = t57 * t69;
t125 = t57 * t71;
t124 = t58 * t57;
t122 = t66 * t57;
t121 = t68 * t57;
t120 = t70 * t57;
t77 = cos(qJ(2,1));
t117 = t77 * t78;
t116 = t55 * t140;
t19 = -t54 * t128 + t56 * t61;
t115 = t19 * t145;
t20 = t54 * t127 - t56 * t73;
t114 = t20 * t144;
t21 = t54 * t126 - t56 * t75;
t113 = t21 * t143;
t64 = legFrame(2,2);
t40 = sin(t64);
t109 = t40 * t121;
t62 = legFrame(4,2);
t38 = sin(t62);
t108 = t38 * t124;
t42 = cos(t62);
t107 = t42 * t124;
t63 = legFrame(3,2);
t39 = sin(t63);
t106 = t39 * t122;
t43 = cos(t63);
t105 = t43 * t122;
t44 = cos(t64);
t104 = t44 * t121;
t103 = (pkin(2) + t139) * t58;
t102 = (pkin(2) + t138) * t66;
t101 = (pkin(2) + t137) * t68;
t33 = t78 * t135;
t23 = t33 + t140;
t24 = pkin(2) * t135 - t129;
t97 = t23 * t61 - t24 * t59;
t96 = t23 * t73 - t24 * t67;
t95 = t23 * t75 - t24 * t69;
t94 = t54 * t103;
t93 = t54 * t102;
t92 = t54 * t101;
t91 = t58 * t141 - t155 * t57;
t90 = t66 * t141 - t154 * t57;
t89 = t68 * t141 - t153 * t57;
t29 = pkin(2) * t71 - t117;
t88 = pkin(3) * t131 - t29 * t57;
t87 = koppelP(1,1);
t86 = koppelP(2,1);
t85 = koppelP(3,1);
t84 = koppelP(4,1);
t83 = koppelP(1,2);
t82 = koppelP(2,2);
t81 = koppelP(3,2);
t80 = koppelP(4,2);
t79 = xP(4);
t65 = legFrame(1,2);
t49 = cos(t79);
t48 = sin(t79);
t45 = cos(t65);
t41 = sin(t65);
t36 = t76 * pkin(3) + pkin(2);
t32 = pkin(2) * t77 + t71 * t78;
t31 = pkin(2) * t75 + t69 * t78;
t30 = pkin(2) * t73 + t67 * t78;
t26 = pkin(2) * t61 + t59 * t78;
t22 = t54 * t125 - t56 * t77;
t18 = pkin(3) * t120 + t55 * t29;
t17 = pkin(3) * t121 + t153 * t55;
t16 = pkin(3) * t122 + t154 * t55;
t15 = pkin(3) * t124 + t155 * t55;
t14 = (t71 * t36 - t117) * t76 * t55 + t36 * t120;
t9 = t56 * t32 + t88 * t54;
t8 = t56 * t31 + t89 * t54;
t7 = t56 * t30 + t90 * t54;
t6 = t56 * t26 + t91 * t54;
t5 = 0.1e1 / (pkin(2) * t120 + t130 * t142 + t18 * t76);
t4 = 0.1e1 / (pkin(2) * t121 + t55 * t110 + t17 * t74);
t3 = 0.1e1 / (pkin(2) * t122 + t55 * t111 + t16 * t72);
t2 = 0.1e1 / (pkin(2) * t124 + t55 * t112 + t15 * t60);
t1 = t54 * t36 * t131 + ((-t36 * t135 + t129) * t71 + (t56 * t36 + t33) * t77) * t76;
t10 = [(-(-t41 * t130 + t22 * t45) * t142 + (t41 * t18 + t9 * t45) * t76 + (t45 * t136 + t57 * t41) * t146) * t5, ((t45 * t130 + t22 * t41) * t142 + (t45 * t18 - t9 * t41) * t76 + (-t41 * t136 + t57 * t45) * t146) * t5, (-(t56 * t125 + t54 * t77) * t142 + (-t32 * t54 + t88 * t56) * t76 + t70 * t116) * t5, (-(t48 * t87 + t49 * t83) * (t1 * t45 + t14 * t41) - (-t48 * t83 + t49 * t87) * (t1 * t41 - t45 * t14)) / t14; (-(-t40 * t132 + t21 * t44) * t143 + (t40 * t17 + t8 * t44) * t74 + (t44 * t136 + t57 * t40) * t147) * t4, ((t44 * t132 + t21 * t40) * t143 + (t44 * t17 - t8 * t40) * t74 + (-t40 * t136 + t57 * t44) * t147) * t4, (-(t56 * t126 + t54 * t75) * t143 + (-t31 * t54 + t89 * t56) * t74 + t68 * t116) * t4, (-(t48 * t86 + t49 * t82) * ((t152 * t40 + t44 * t92) * t55 - t44 * t113 + (pkin(3) * t109 + t95 * t44) * t74 + pkin(2) * t109) - (-t48 * t82 + t49 * t86) * ((-t152 * t44 + t40 * t92) * t55 - t40 * t113 + (-pkin(3) * t104 + t95 * t40) * t74 - pkin(2) * t104)) / ((t69 * t137 + t153) * t74 * t55 + t57 * t101); (-(-t39 * t133 + t20 * t43) * t144 + (t39 * t16 + t7 * t43) * t72 + (t43 * t136 + t57 * t39) * t148) * t3, ((t43 * t133 + t20 * t39) * t144 + (t43 * t16 - t7 * t39) * t72 + (-t39 * t136 + t57 * t43) * t148) * t3, (-(t56 * t127 + t54 * t73) * t144 + (-t30 * t54 + t90 * t56) * t72 + t66 * t116) * t3, (-(t48 * t85 + t49 * t81) * ((t151 * t39 + t43 * t93) * t55 - t43 * t114 + (pkin(3) * t106 + t96 * t43) * t72 + pkin(2) * t106) - (-t48 * t81 + t49 * t85) * ((-t151 * t43 + t39 * t93) * t55 - t39 * t114 + (-pkin(3) * t105 + t96 * t39) * t72 - pkin(2) * t105)) / ((t67 * t138 + t154) * t72 * t55 + t57 * t102); ((t38 * t134 + t19 * t42) * t145 + (t38 * t15 + t6 * t42) * t60 + (t42 * t136 + t57 * t38) * t149) * t2, (-(-t42 * t134 + t19 * t38) * t145 + (t42 * t15 - t6 * t38) * t60 + (-t38 * t136 + t57 * t42) * t149) * t2, (-(t56 * t128 + t54 * t61) * t145 + (-t26 * t54 + t91 * t56) * t60 + t58 * t116) * t2, (-(t48 * t84 + t49 * t80) * ((t150 * t38 + t42 * t94) * t55 + t42 * t115 + (pkin(3) * t108 + t97 * t42) * t60 + pkin(2) * t108) - ((-t150 * t42 + t38 * t94) * t55 + t38 * t115 + (-pkin(3) * t107 + t97 * t38) * t60 - pkin(2) * t107) * (-t48 * t80 + t49 * t84)) / ((t59 * t139 + t155) * t60 * t55 + t57 * t103);];
Jinv  = t10;
