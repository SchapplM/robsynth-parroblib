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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d3,d4,theta1,theta2]';
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
% Datum: 2019-05-16 22:39
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P6PPRRRR3V2A3_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(10,1),zeros(6,3),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPRRRR3V2A3_Jinv: qJ has to be [3x6] (double)');
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPRRRR3V2A3_Jinv: xP has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPRRRR3V2A3_Jinv: pkin has to be [10x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPRRRR3V2A3_Jinv: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPRRRR3V2A3_Jinv: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 22:37:45
% EndTime: 2019-05-16 22:37:47
% DurationCPUTime: 1.63s
% Computational Cost: add. (642->222), mult. (1566->500), div. (36->6), fcn. (1626->60), ass. (0->231)
t139 = legFrame(6,2);
t103 = cos(t139);
t151 = xP(6);
t109 = sin(t151);
t152 = xP(5);
t110 = sin(t152);
t112 = cos(t151);
t153 = xP(4);
t114 = cos(t153);
t166 = koppelP(6,1);
t208 = t114 * t166;
t160 = koppelP(6,2);
t214 = t114 * t160;
t113 = cos(t152);
t154 = koppelP(6,3);
t220 = t113 * t154;
t111 = sin(t153);
t226 = t111 * t166;
t232 = t111 * t160;
t189 = t109 * t160 - t112 * t166;
t25 = t110 * t154 - t189 * t113;
t127 = legFrame(6,1);
t91 = cos(t127);
t97 = sin(t139);
t265 = t91 * t97;
t281 = t103 * (t109 * (t110 * t214 + t226) - (t110 * t208 - t232) * t112 + t114 * t220) + t25 * t265;
t79 = sin(t127);
t268 = t79 * t97;
t280 = t103 * (-t109 * (-t110 * t232 + t208) + t111 * t220 - (t110 * t226 + t214) * t112) + t25 * t268;
t140 = legFrame(5,2);
t104 = cos(t140);
t167 = koppelP(5,1);
t207 = t114 * t167;
t161 = koppelP(5,2);
t213 = t114 * t161;
t155 = koppelP(5,3);
t219 = t113 * t155;
t225 = t111 * t167;
t231 = t111 * t161;
t188 = t109 * t161 - t112 * t167;
t26 = t110 * t155 - t188 * t113;
t128 = legFrame(5,1);
t92 = cos(t128);
t98 = sin(t140);
t264 = t92 * t98;
t279 = t104 * (t109 * (t110 * t213 + t225) - (t110 * t207 - t231) * t112 + t114 * t219) + t26 * t264;
t80 = sin(t128);
t267 = t80 * t98;
t278 = t104 * (-t109 * (-t110 * t231 + t207) + t111 * t219 - (t110 * t225 + t213) * t112) + t26 * t267;
t141 = legFrame(4,2);
t105 = cos(t141);
t168 = koppelP(4,1);
t206 = t114 * t168;
t162 = koppelP(4,2);
t212 = t114 * t162;
t156 = koppelP(4,3);
t218 = t113 * t156;
t224 = t111 * t168;
t230 = t111 * t162;
t129 = legFrame(4,1);
t93 = cos(t129);
t99 = sin(t141);
t263 = t93 * t99;
t187 = t109 * t162 - t112 * t168;
t27 = t110 * t156 - t187 * t113;
t277 = t105 * (t109 * (t110 * t212 + t224) - (t110 * t206 - t230) * t112 + t114 * t218) + t27 * t263;
t81 = sin(t129);
t266 = t81 * t99;
t276 = t105 * (-t109 * (-t110 * t230 + t206) + t111 * t218 - (t110 * t224 + t212) * t112) + t27 * t266;
t142 = legFrame(3,2);
t106 = cos(t142);
t169 = koppelP(3,1);
t205 = t114 * t169;
t163 = koppelP(3,2);
t211 = t114 * t163;
t157 = koppelP(3,3);
t217 = t113 * t157;
t223 = t111 * t169;
t229 = t111 * t163;
t100 = sin(t142);
t130 = legFrame(3,1);
t94 = cos(t130);
t241 = t94 * t100;
t186 = t109 * t163 - t112 * t169;
t28 = t110 * t157 - t186 * t113;
t275 = t106 * (t109 * (t110 * t211 + t223) - (t110 * t205 - t229) * t112 + t114 * t217) + t28 * t241;
t82 = sin(t130);
t244 = t82 * t100;
t274 = t106 * (-t109 * (-t110 * t229 + t205) + t111 * t217 - (t110 * t223 + t211) * t112) + t28 * t244;
t143 = legFrame(2,2);
t107 = cos(t143);
t170 = koppelP(2,1);
t204 = t114 * t170;
t164 = koppelP(2,2);
t210 = t114 * t164;
t158 = koppelP(2,3);
t216 = t113 * t158;
t222 = t111 * t170;
t228 = t111 * t164;
t101 = sin(t143);
t131 = legFrame(2,1);
t95 = cos(t131);
t240 = t95 * t101;
t185 = t109 * t164 - t112 * t170;
t29 = t110 * t158 - t185 * t113;
t273 = t107 * (t109 * (t110 * t210 + t222) - (t110 * t204 - t228) * t112 + t114 * t216) + t29 * t240;
t83 = sin(t131);
t243 = t83 * t101;
t272 = t107 * (-t109 * (-t110 * t228 + t204) + t111 * t216 - (t110 * t222 + t210) * t112) + t29 * t243;
t144 = legFrame(1,2);
t108 = cos(t144);
t171 = koppelP(1,1);
t203 = t114 * t171;
t165 = koppelP(1,2);
t209 = t114 * t165;
t159 = koppelP(1,3);
t215 = t113 * t159;
t221 = t111 * t171;
t227 = t111 * t165;
t102 = sin(t144);
t132 = legFrame(1,1);
t96 = cos(t132);
t239 = t96 * t102;
t184 = t109 * t165 - t112 * t171;
t30 = t110 * t159 - t184 * t113;
t271 = (t109 * (t110 * t209 + t221) - (t110 * t203 - t227) * t112 + t114 * t215) * t108 + t30 * t239;
t84 = sin(t132);
t242 = t84 * t102;
t270 = t108 * (-t109 * (-t110 * t227 + t203) + t111 * t215 - (t110 * t221 + t209) * t112) + t30 * t242;
t269 = pkin(8) * sin(pkin(6));
t116 = sin(pkin(9));
t262 = t116 * t25;
t261 = t116 * t26;
t260 = t116 * t27;
t259 = t116 * t28;
t258 = t116 * t29;
t257 = t116 * t30;
t256 = t116 * t97;
t255 = t116 * t98;
t254 = t116 * t99;
t119 = cos(pkin(9));
t253 = t119 * t97;
t252 = t119 * t98;
t251 = t119 * t99;
t250 = t25 * t119;
t249 = t26 * t119;
t248 = t27 * t119;
t247 = t28 * t119;
t246 = t29 * t119;
t245 = t30 * t119;
t238 = t100 * t116;
t237 = t100 * t119;
t236 = t101 * t116;
t235 = t101 * t119;
t234 = t102 * t116;
t233 = t102 * t119;
t202 = sin(pkin(10)) * cos(pkin(5));
t150 = cos(qJ(3,1));
t149 = cos(qJ(3,2));
t148 = cos(qJ(3,3));
t147 = sin(qJ(3,1));
t146 = sin(qJ(3,2));
t145 = sin(qJ(3,3));
t138 = cos(qJ(3,4));
t137 = cos(qJ(3,5));
t136 = cos(qJ(3,6));
t135 = sin(qJ(3,4));
t134 = sin(qJ(3,5));
t133 = sin(qJ(3,6));
t126 = legFrame(1,3);
t125 = legFrame(2,3);
t124 = legFrame(3,3);
t123 = legFrame(4,3);
t122 = legFrame(5,3);
t121 = legFrame(6,3);
t118 = cos(pkin(10));
t90 = cos(t126);
t89 = cos(t125);
t88 = cos(t124);
t87 = cos(t123);
t86 = cos(t122);
t85 = cos(t121);
t78 = sin(t126);
t77 = sin(t125);
t76 = sin(t124);
t75 = sin(t123);
t74 = sin(t122);
t73 = sin(t121);
t72 = t109 * t171 + t112 * t165;
t71 = t109 * t170 + t112 * t164;
t70 = t109 * t169 + t112 * t163;
t69 = t109 * t168 + t112 * t162;
t68 = t109 * t167 + t112 * t161;
t67 = t109 * t166 + t112 * t160;
t66 = -t78 * t116 + t119 * t90;
t65 = -t77 * t116 + t119 * t89;
t64 = -t76 * t116 + t119 * t88;
t63 = -t75 * t116 + t119 * t87;
t62 = -t74 * t116 + t119 * t86;
t61 = -t73 * t116 + t119 * t85;
t60 = t116 * t90 + t119 * t78;
t59 = t116 * t89 + t119 * t77;
t58 = t116 * t88 + t119 * t76;
t57 = t116 * t87 + t119 * t75;
t56 = t116 * t86 + t119 * t74;
t55 = t116 * t85 + t119 * t73;
t24 = t184 * t110 + t215;
t23 = t185 * t110 + t216;
t22 = t186 * t110 + t217;
t21 = t187 * t110 + t218;
t20 = t188 * t110 + t219;
t19 = t189 * t110 + t220;
t18 = 0.1e1 / ((pkin(3) * t150 + t147 * t269) * t202 + t118 * (pkin(3) * t147 - t150 * t269));
t17 = 0.1e1 / ((pkin(3) * t149 + t146 * t269) * t202 + t118 * (pkin(3) * t146 - t149 * t269));
t16 = 0.1e1 / ((pkin(3) * t148 + t145 * t269) * t202 + t118 * (pkin(3) * t145 - t148 * t269));
t15 = 0.1e1 / ((pkin(3) * t138 + t135 * t269) * t202 + t118 * (pkin(3) * t135 - t138 * t269));
t14 = 0.1e1 / ((pkin(3) * t137 + t134 * t269) * t202 + t118 * (pkin(3) * t134 - t137 * t269));
t13 = 0.1e1 / ((pkin(3) * t136 + t133 * t269) * t202 + t118 * (pkin(3) * t133 - t136 * t269));
t12 = t111 * t24 - t72 * t114;
t11 = t111 * t23 - t71 * t114;
t10 = t111 * t22 - t70 * t114;
t9 = t111 * t21 - t69 * t114;
t8 = t111 * t20 - t68 * t114;
t7 = t111 * t19 - t67 * t114;
t6 = t111 * t72 + t24 * t114;
t5 = t111 * t71 + t23 * t114;
t4 = t111 * t70 + t22 * t114;
t3 = t111 * t69 + t21 * t114;
t2 = t111 * t68 + t20 * t114;
t1 = t111 * t67 + t19 * t114;
t31 = [-t66 * t108 * t18, (-t66 * t242 - t60 * t96) * t18, (t66 * t239 - t84 * t60) * t18, (((t116 * t6 - t12 * t233) * t90 + t78 * (t119 * t6 + t12 * t234)) * t96 + ((t116 * t12 + t6 * t233) * t90 + t78 * (t119 * t12 - t6 * t234)) * t84) * t18, ((-t271 * t119 + t84 * t257) * t90 + (t116 * t271 + t84 * t245) * t78) * t18, ((-t270 * t119 - t96 * t257) * t90 - (-t116 * t270 + t96 * t245) * t78) * t18; -t65 * t107 * t17, (-t65 * t243 - t59 * t95) * t17, (t65 * t240 - t83 * t59) * t17, (((-t11 * t235 + t116 * t5) * t89 + t77 * (t11 * t236 + t119 * t5)) * t95 + ((t116 * t11 + t5 * t235) * t89 + t77 * (t119 * t11 - t5 * t236)) * t83) * t17, ((-t273 * t119 + t83 * t258) * t89 + (t116 * t273 + t83 * t246) * t77) * t17, ((-t272 * t119 - t95 * t258) * t89 - (-t116 * t272 + t95 * t246) * t77) * t17; -t64 * t106 * t16, (-t64 * t244 - t58 * t94) * t16, (t64 * t241 - t82 * t58) * t16, (((-t10 * t237 + t116 * t4) * t88 + t76 * (t10 * t238 + t119 * t4)) * t94 + ((t116 * t10 + t4 * t237) * t88 + t76 * (t119 * t10 - t4 * t238)) * t82) * t16, ((-t275 * t119 + t82 * t259) * t88 + (t116 * t275 + t82 * t247) * t76) * t16, ((-t274 * t119 - t94 * t259) * t88 - (-t116 * t274 + t94 * t247) * t76) * t16; -t63 * t105 * t15, (-t63 * t266 - t57 * t93) * t15, (t63 * t263 - t81 * t57) * t15, (((t116 * t3 - t9 * t251) * t87 + t75 * (t119 * t3 + t9 * t254)) * t93 + ((t116 * t9 + t3 * t251) * t87 + t75 * (t119 * t9 - t3 * t254)) * t81) * t15, ((-t277 * t119 + t81 * t260) * t87 + (t116 * t277 + t81 * t248) * t75) * t15, ((-t276 * t119 - t93 * t260) * t87 - (-t116 * t276 + t93 * t248) * t75) * t15; -t62 * t104 * t14, (-t62 * t267 - t56 * t92) * t14, (t62 * t264 - t80 * t56) * t14, (((t116 * t2 - t8 * t252) * t86 + t74 * (t119 * t2 + t8 * t255)) * t92 + ((t116 * t8 + t2 * t252) * t86 + t74 * (t119 * t8 - t2 * t255)) * t80) * t14, ((-t279 * t119 + t80 * t261) * t86 + (t116 * t279 + t80 * t249) * t74) * t14, ((-t278 * t119 - t92 * t261) * t86 - (-t116 * t278 + t92 * t249) * t74) * t14; -t61 * t103 * t13, (-t61 * t268 - t55 * t91) * t13, (t61 * t265 - t79 * t55) * t13, (((t116 * t1 - t7 * t253) * t85 + t73 * (t119 * t1 + t7 * t256)) * t91 + ((t1 * t253 + t116 * t7) * t85 + t73 * (-t1 * t256 + t119 * t7)) * t79) * t13, ((-t281 * t119 + t79 * t262) * t85 + (t116 * t281 + t79 * t250) * t73) * t13, ((-t280 * t119 - t91 * t262) * t85 - (-t116 * t280 + t91 * t250) * t73) * t13;];
Jinv  = t31;
