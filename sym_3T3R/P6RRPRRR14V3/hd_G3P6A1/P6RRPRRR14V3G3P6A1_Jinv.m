% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-09-28 23:14
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V3G3P6A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(1,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G3P6A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G3P6A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G3P6A1_Jinv: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G3P6A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G3P6A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-09-28 23:14:29
% EndTime: 2020-09-28 23:14:30
% DurationCPUTime: 0.38s
% Computational Cost: add. (180->102), mult. (456->186), div. (0->0), fcn. (510->42), ass. (0->130)
t81 = xP(4);
t45 = sin(t81);
t88 = koppelP(6,2);
t147 = t45 * t88;
t89 = koppelP(5,2);
t146 = t45 * t89;
t90 = koppelP(4,2);
t145 = t45 * t90;
t91 = koppelP(3,2);
t144 = t45 * t91;
t92 = koppelP(2,2);
t143 = t45 * t92;
t93 = koppelP(1,2);
t142 = t45 * t93;
t94 = koppelP(6,1);
t141 = t45 * t94;
t95 = koppelP(5,1);
t140 = t45 * t95;
t96 = koppelP(4,1);
t139 = t45 * t96;
t97 = koppelP(3,1);
t138 = t45 * t97;
t98 = koppelP(2,1);
t137 = t45 * t98;
t99 = koppelP(1,1);
t136 = t45 * t99;
t80 = xP(5);
t47 = cos(t80);
t82 = koppelP(6,3);
t135 = t47 * t82;
t83 = koppelP(5,3);
t134 = t47 * t83;
t84 = koppelP(4,3);
t133 = t47 * t84;
t85 = koppelP(3,3);
t132 = t47 * t85;
t86 = koppelP(2,3);
t131 = t47 * t86;
t87 = koppelP(1,3);
t130 = t47 * t87;
t48 = cos(t81);
t129 = t48 * t88;
t128 = t48 * t89;
t127 = t48 * t90;
t126 = t48 * t91;
t125 = t48 * t92;
t124 = t48 * t93;
t123 = t48 * t94;
t122 = t48 * t95;
t121 = t48 * t96;
t120 = t48 * t97;
t119 = t48 * t98;
t118 = t48 * t99;
t49 = sin(qJ(2,6));
t117 = t49 * cos(qJ(1,6));
t116 = sin(qJ(1,6)) * t49;
t51 = sin(qJ(2,5));
t115 = t51 * cos(qJ(1,5));
t114 = sin(qJ(1,5)) * t51;
t53 = sin(qJ(2,4));
t113 = t53 * cos(qJ(1,4));
t112 = sin(qJ(1,4)) * t53;
t67 = sin(qJ(2,3));
t111 = t67 * cos(qJ(1,3));
t110 = sin(qJ(1,3)) * t67;
t69 = sin(qJ(2,2));
t109 = t69 * cos(qJ(1,2));
t108 = sin(qJ(1,2)) * t69;
t71 = sin(qJ(2,1));
t107 = t71 * cos(qJ(1,1));
t106 = sin(qJ(1,1)) * t71;
t79 = xP(6);
t43 = sin(t79);
t44 = sin(t80);
t46 = cos(t79);
t105 = (-t43 * t88 + t46 * t94) * t47 + t44 * t82;
t104 = (-t43 * t89 + t46 * t95) * t47 + t44 * t83;
t103 = (-t43 * t90 + t46 * t96) * t47 + t44 * t84;
t102 = (-t43 * t91 + t46 * t97) * t47 + t44 * t85;
t101 = (-t43 * t92 + t46 * t98) * t47 + t44 * t86;
t100 = (-t43 * t93 + t46 * t99) * t47 + t44 * t87;
t77 = cos(qJ(2,1));
t75 = cos(qJ(2,2));
t73 = cos(qJ(2,3));
t66 = legFrame(1,2);
t65 = legFrame(2,2);
t64 = legFrame(3,2);
t63 = legFrame(4,2);
t62 = legFrame(5,2);
t61 = legFrame(6,2);
t59 = cos(qJ(2,4));
t57 = cos(qJ(2,5));
t55 = cos(qJ(2,6));
t42 = cos(t66);
t41 = cos(t65);
t40 = cos(t64);
t39 = cos(t63);
t38 = cos(t62);
t37 = cos(t61);
t36 = sin(t66);
t35 = sin(t65);
t34 = sin(t64);
t33 = sin(t63);
t32 = sin(t62);
t31 = sin(t61);
t24 = t42 * t107 - t36 * t77;
t23 = t41 * t109 - t35 * t75;
t22 = t40 * t111 - t34 * t73;
t21 = t36 * t107 + t42 * t77;
t20 = t35 * t109 + t41 * t75;
t19 = t34 * t111 + t40 * t73;
t18 = t39 * t113 - t33 * t59;
t17 = t38 * t115 - t32 * t57;
t16 = t37 * t117 - t31 * t55;
t15 = t33 * t113 + t39 * t59;
t14 = t32 * t115 + t38 * t57;
t13 = t31 * t117 + t37 * t55;
t12 = (t44 * t121 - t145) * t46 + (-t44 * t127 - t139) * t43 - t48 * t133;
t11 = (t44 * t122 - t146) * t46 + (-t44 * t128 - t140) * t43 - t48 * t134;
t10 = (t44 * t123 - t147) * t46 + (-t44 * t129 - t141) * t43 - t48 * t135;
t9 = (-t44 * t124 - t136) * t43 + (t44 * t118 - t142) * t46 - t48 * t130;
t8 = (-t44 * t142 + t118) * t43 + (t44 * t136 + t124) * t46 - t45 * t130;
t7 = (-t44 * t125 - t137) * t43 + (t44 * t119 - t143) * t46 - t48 * t131;
t6 = (-t44 * t126 - t138) * t43 + (t44 * t120 - t144) * t46 - t48 * t132;
t5 = -t45 * t131 + (t44 * t137 + t125) * t46 + t43 * (-t44 * t143 + t119);
t4 = -t45 * t132 + (t44 * t138 + t126) * t46 + t43 * (-t44 * t144 + t120);
t3 = -t45 * t133 + (t44 * t139 + t127) * t46 + t43 * (-t44 * t145 + t121);
t2 = -t45 * t134 + (t44 * t140 + t128) * t46 + t43 * (-t44 * t146 + t122);
t1 = -t45 * t135 + (t44 * t141 + t129) * t46 + t43 * (-t44 * t147 + t123);
t25 = [t24, -t21, -t106, -t8 * t106 - t21 * t9, t100 * t106 - t24 * t9, -t21 * t100 - t24 * t8; t23, -t20, -t108, -t5 * t108 - t20 * t7, t101 * t108 - t23 * t7, -t20 * t101 - t23 * t5; t22, -t19, -t110, -t4 * t110 - t19 * t6, t102 * t110 - t22 * t6, -t19 * t102 - t22 * t4; t18, -t15, -t112, -t3 * t112 - t12 * t15, t103 * t112 - t12 * t18, -t15 * t103 - t3 * t18; t17, -t14, -t114, -t11 * t14 - t2 * t114, t104 * t114 - t11 * t17, -t14 * t104 - t2 * t17; t16, -t13, -t116, -t1 * t116 - t10 * t13, -t10 * t16 + t105 * t116, -t1 * t16 - t13 * t105;];
Jinv  = t25;
