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
% Datum: 2020-08-23 02:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6RRPRRR14V3G1P4A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(1,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1P4A1_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1P4A1_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1P4A1_Jinv: pkin has to be [1x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1P4A1_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1P4A1_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-23 02:42:49
% EndTime: 2020-08-23 02:42:49
% DurationCPUTime: 0.42s
% Computational Cost: add. (216->114), mult. (498->217), div. (0->0), fcn. (540->42), ass. (0->125)
t73 = legFrame(6,3);
t55 = sin(t73);
t61 = cos(t73);
t79 = sin(qJ(2,6));
t80 = sin(qJ(1,6));
t86 = cos(qJ(1,6));
t154 = (t55 * t80 - t86 * t61) * t79;
t74 = legFrame(5,3);
t56 = sin(t74);
t62 = cos(t74);
t81 = sin(qJ(2,5));
t82 = sin(qJ(1,5));
t88 = cos(qJ(1,5));
t153 = (t56 * t82 - t88 * t62) * t81;
t75 = legFrame(4,3);
t57 = sin(t75);
t63 = cos(t75);
t83 = sin(qJ(2,4));
t84 = sin(qJ(1,4));
t90 = cos(qJ(1,4));
t152 = (t57 * t84 - t90 * t63) * t83;
t76 = legFrame(3,3);
t58 = sin(t76);
t64 = cos(t76);
t91 = sin(qJ(2,3));
t92 = sin(qJ(1,3));
t98 = cos(qJ(1,3));
t151 = (t58 * t92 - t98 * t64) * t91;
t100 = cos(qJ(1,2));
t77 = legFrame(2,3);
t59 = sin(t77);
t65 = cos(t77);
t93 = sin(qJ(2,2));
t94 = sin(qJ(1,2));
t150 = (-t100 * t65 + t59 * t94) * t93;
t102 = cos(qJ(1,1));
t78 = legFrame(1,3);
t60 = sin(t78);
t66 = cos(t78);
t95 = sin(qJ(2,1));
t96 = sin(qJ(1,1));
t149 = (-t102 * t66 + t60 * t96) * t95;
t104 = xP(5);
t68 = sin(t104);
t105 = xP(4);
t69 = sin(t105);
t148 = t68 * t69;
t147 = t79 * (t55 * t86 + t80 * t61);
t146 = t81 * (t56 * t88 + t82 * t62);
t145 = t83 * (t57 * t90 + t84 * t63);
t144 = t91 * (t58 * t98 + t92 * t64);
t143 = t93 * (t59 * t100 + t94 * t65);
t142 = t95 * (t60 * t102 + t96 * t66);
t106 = koppelP(6,3);
t71 = cos(t104);
t141 = t106 * t71;
t107 = koppelP(5,3);
t140 = t107 * t71;
t108 = koppelP(4,3);
t139 = t108 * t71;
t109 = koppelP(3,3);
t138 = t109 * t71;
t110 = koppelP(2,3);
t137 = t110 * t71;
t111 = koppelP(1,3);
t136 = t111 * t71;
t112 = koppelP(6,2);
t118 = koppelP(6,1);
t103 = xP(6);
t67 = sin(t103);
t70 = cos(t103);
t135 = t112 * t67 - t118 * t70;
t113 = koppelP(5,2);
t119 = koppelP(5,1);
t134 = t113 * t67 - t119 * t70;
t114 = koppelP(4,2);
t120 = koppelP(4,1);
t133 = t114 * t67 - t120 * t70;
t115 = koppelP(3,2);
t121 = koppelP(3,1);
t132 = t115 * t67 - t121 * t70;
t116 = koppelP(2,2);
t122 = koppelP(2,1);
t131 = t116 * t67 - t122 * t70;
t117 = koppelP(1,2);
t123 = koppelP(1,1);
t130 = t117 * t67 - t123 * t70;
t13 = t68 * t106 - t135 * t71;
t14 = t68 * t107 - t134 * t71;
t15 = t68 * t108 - t133 * t71;
t16 = t68 * t109 - t132 * t71;
t17 = t68 * t110 - t131 * t71;
t18 = t68 * t111 - t130 * t71;
t72 = cos(t105);
t129 = t69 * t141 - t67 * (-t112 * t148 + t118 * t72) - (t112 * t72 + t118 * t148) * t70;
t128 = t69 * t140 - t67 * (-t113 * t148 + t119 * t72) - (t113 * t72 + t119 * t148) * t70;
t127 = t69 * t139 - t67 * (-t114 * t148 + t120 * t72) - (t114 * t72 + t120 * t148) * t70;
t126 = t69 * t138 - t67 * (-t115 * t148 + t121 * t72) - (t115 * t72 + t121 * t148) * t70;
t125 = t69 * t137 - t67 * (-t116 * t148 + t122 * t72) - (t116 * t72 + t122 * t148) * t70;
t124 = t69 * t136 - t67 * (-t117 * t148 + t123 * t72) - (t117 * t72 + t123 * t148) * t70;
t101 = cos(qJ(2,1));
t99 = cos(qJ(2,2));
t97 = cos(qJ(2,3));
t89 = cos(qJ(2,4));
t87 = cos(qJ(2,5));
t85 = cos(qJ(2,6));
t48 = t117 * t70 + t123 * t67;
t47 = t116 * t70 + t122 * t67;
t46 = t115 * t70 + t121 * t67;
t45 = t114 * t70 + t120 * t67;
t44 = t113 * t70 + t119 * t67;
t43 = t112 * t70 + t118 * t67;
t12 = t130 * t68 + t136;
t11 = t131 * t68 + t137;
t10 = t132 * t68 + t138;
t9 = t133 * t68 + t139;
t8 = t134 * t68 + t140;
t7 = t135 * t68 + t141;
t6 = t12 * t72 + t69 * t48;
t5 = t11 * t72 + t69 * t47;
t4 = t10 * t72 + t69 * t46;
t3 = t69 * t45 + t9 * t72;
t2 = t69 * t44 + t8 * t72;
t1 = t69 * t43 + t7 * t72;
t19 = [-t149, t142, -t101, -t6 * t142 - t101 * (-t69 * t12 + t72 * t48), t101 * t18 - t6 * t149, t95 * ((t124 * t102 + t96 * t18) * t66 + (t18 * t102 - t96 * t124) * t60); -t150, t143, -t99, -t5 * t143 - t99 * (-t69 * t11 + t72 * t47), -t5 * t150 + t99 * t17, t93 * ((t125 * t100 + t94 * t17) * t65 + (t17 * t100 - t94 * t125) * t59); -t151, t144, -t97, -t4 * t144 - t97 * (-t69 * t10 + t72 * t46), -t4 * t151 + t97 * t16, t91 * ((t126 * t98 + t92 * t16) * t64 + (-t92 * t126 + t16 * t98) * t58); -t152, t145, -t89, -t3 * t145 - t89 * (t72 * t45 - t69 * t9), t89 * t15 - t3 * t152, t83 * ((t127 * t90 + t84 * t15) * t63 + (-t84 * t127 + t15 * t90) * t57); -t153, t146, -t87, -t2 * t146 - t87 * (t72 * t44 - t69 * t8), t87 * t14 - t2 * t153, t81 * ((t128 * t88 + t82 * t14) * t62 + (-t82 * t128 + t14 * t88) * t56); -t154, t147, -t85, -t1 * t147 - t85 * (t72 * t43 - t69 * t7), -t1 * t154 + t85 * t13, t79 * ((t129 * t86 + t80 * t13) * t61 + (-t80 * t129 + t13 * t86) * t55);];
Jinv  = t19;
