% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:11
% EndTime: 2020-08-06 19:59:15
% DurationCPUTime: 4.56s
% Computational Cost: add. (16977->426), mult. (28131->723), div. (3315->12), fcn. (23382->35), ass. (0->277)
t172 = pkin(4) + qJ(3,3);
t164 = 0.1e1 / t172;
t196 = xDP(3);
t197 = xDP(2);
t198 = xDP(1);
t177 = legFrame(3,2);
t151 = sin(t177);
t154 = cos(t177);
t184 = sin(qJ(1,3));
t190 = cos(qJ(1,3));
t171 = cos(pkin(5));
t289 = pkin(2) * t171;
t131 = pkin(1) + t289;
t189 = cos(qJ(2,3));
t170 = sin(pkin(5));
t183 = sin(qJ(2,3));
t258 = t170 * t183;
t215 = pkin(2) * t258 - t131 * t189;
t83 = -t172 * t190 - t184 * t215;
t290 = pkin(2) * t170;
t96 = t131 * t183 + t189 * t290;
t55 = t151 * t96 + t154 * t83;
t56 = -t151 * t83 + t154 * t96;
t157 = t189 * pkin(1);
t161 = qJ(2,3) + pkin(5);
t318 = t157 + pkin(2) * cos(t161);
t87 = t172 * t184 + t190 * t318;
t40 = (t196 * t87 + t197 * t56 + t198 * t55) * t164;
t328 = 0.2e1 * t40;
t173 = pkin(4) + qJ(3,2);
t165 = 0.1e1 / t173;
t178 = legFrame(2,2);
t152 = sin(t178);
t155 = cos(t178);
t186 = sin(qJ(1,2));
t192 = cos(qJ(1,2));
t191 = cos(qJ(2,2));
t185 = sin(qJ(2,2));
t257 = t170 * t185;
t214 = pkin(2) * t257 - t131 * t191;
t84 = -t173 * t192 - t186 * t214;
t97 = t131 * t185 + t191 * t290;
t57 = t152 * t97 + t155 * t84;
t58 = -t152 * t84 + t155 * t97;
t158 = t191 * pkin(1);
t162 = qJ(2,2) + pkin(5);
t317 = t158 + pkin(2) * cos(t162);
t88 = t173 * t186 + t192 * t317;
t41 = (t196 * t88 + t197 * t58 + t198 * t57) * t165;
t327 = 0.2e1 * t41;
t174 = pkin(4) + qJ(3,1);
t166 = 0.1e1 / t174;
t179 = legFrame(1,2);
t153 = sin(t179);
t156 = cos(t179);
t188 = sin(qJ(1,1));
t194 = cos(qJ(1,1));
t193 = cos(qJ(2,1));
t187 = sin(qJ(2,1));
t256 = t170 * t187;
t213 = pkin(2) * t256 - t131 * t193;
t85 = -t174 * t194 - t188 * t213;
t98 = t131 * t187 + t193 * t290;
t59 = t153 * t98 + t156 * t85;
t60 = -t153 * t85 + t156 * t98;
t159 = t193 * pkin(1);
t163 = qJ(2,1) + pkin(5);
t316 = t159 + pkin(2) * cos(t163);
t89 = t174 * t188 + t194 * t316;
t42 = (t196 * t89 + t197 * t60 + t198 * t59) * t166;
t326 = 0.2e1 * t42;
t145 = t170 * mrSges(3,1);
t319 = t171 * mrSges(3,2) + t145;
t253 = 0.2e1 * pkin(1);
t305 = 0.2e1 * t171;
t325 = t305 / 0.2e1;
t323 = mrSges(2,3) - mrSges(1,2);
t203 = pkin(2) ^ 2;
t204 = pkin(1) ^ 2;
t254 = -t203 - t204;
t121 = t289 * t253 - t254;
t109 = 0.1e1 / t318;
t73 = (t151 * t198 + t154 * t197) * t109;
t70 = t73 ^ 2;
t322 = t70 * t121;
t110 = 0.1e1 / t317;
t74 = (t152 * t198 + t155 * t197) * t110;
t71 = t74 ^ 2;
t321 = t71 * t121;
t111 = 0.1e1 / t316;
t75 = (t153 * t198 + t156 * t197) * t111;
t72 = t75 ^ 2;
t320 = t72 * t121;
t315 = pkin(1) * t145 - Ifges(2,4) + Ifges(3,4);
t141 = mrSges(3,2) * qJ(3,1) - Ifges(3,6);
t144 = mrSges(3,1) * qJ(3,1) - Ifges(3,5);
t300 = m(3) * qJ(3,1);
t150 = mrSges(3,3) + t300;
t314 = -pkin(1) * t150 + t141 * t170 - t144 * t171 + Ifges(2,5);
t140 = mrSges(3,2) * qJ(3,2) - Ifges(3,6);
t143 = mrSges(3,1) * qJ(3,2) - Ifges(3,5);
t299 = m(3) * qJ(3,2);
t149 = mrSges(3,3) + t299;
t313 = -pkin(1) * t149 + t140 * t170 - t143 * t171 + Ifges(2,5);
t139 = mrSges(3,2) * qJ(3,3) - Ifges(3,6);
t142 = mrSges(3,1) * qJ(3,3) - Ifges(3,5);
t298 = m(3) * qJ(3,3);
t148 = mrSges(3,3) + t298;
t312 = -pkin(1) * t148 + t139 * t170 - t142 * t171 + Ifges(2,5);
t311 = 0.2e1 * t319;
t175 = Ifges(3,2) - Ifges(3,1);
t279 = mrSges(3,2) * t170;
t297 = m(3) * t204;
t101 = -0.2e1 * pkin(1) * t279 - Ifges(2,1) + Ifges(2,2) - t175 + t297;
t160 = t171 ^ 2;
t255 = t175 * t160;
t277 = Ifges(3,4) * t170;
t301 = pkin(1) * mrSges(3,1);
t76 = (0.4e1 * t277 + 0.2e1 * t301) * t171 + t101 + 0.2e1 * t255;
t309 = 0.2e1 * mrSges(3,3);
t246 = t184 * t290;
t261 = t131 * t184;
t61 = (-t151 * t261 + t154 * t290) * t189 + t183 * (t131 * t154 + t151 * t246);
t64 = (t151 * t290 + t154 * t261) * t189 + (t131 * t151 - t154 * t246) * t183;
t93 = 0.1e1 / t215;
t28 = (t190 * t196 - (t197 * t61 + t198 * t64) * t93) * t164;
t308 = 0.2e1 * t28;
t245 = t186 * t290;
t260 = t131 * t186;
t62 = (-t152 * t260 + t155 * t290) * t191 + t185 * (t131 * t155 + t152 * t245);
t65 = (t152 * t290 + t155 * t260) * t191 + (t131 * t152 - t155 * t245) * t185;
t94 = 0.1e1 / t214;
t29 = (t192 * t196 - (t197 * t62 + t198 * t65) * t94) * t165;
t307 = 0.2e1 * t29;
t244 = t188 * t290;
t259 = t131 * t188;
t63 = (-t153 * t259 + t156 * t290) * t193 + t187 * (t131 * t156 + t153 * t244);
t66 = (t153 * t290 + t156 * t259) * t193 + (t131 * t153 - t156 * t244) * t187;
t95 = 0.1e1 / t213;
t30 = (t194 * t196 - (t197 * t63 + t198 * t66) * t95) * t166;
t306 = 0.2e1 * t30;
t304 = -0.2e1 * t183;
t303 = -0.2e1 * t185;
t302 = -0.2e1 * t187;
t288 = t61 * t93;
t287 = t62 * t94;
t286 = t63 * t95;
t285 = t64 * t93;
t284 = t65 * t94;
t283 = t66 * t95;
t125 = m(3) * pkin(1) - t279;
t147 = mrSges(3,1) * t171;
t102 = -t125 - t147;
t77 = t102 * t189 + t183 * t319;
t282 = t77 * t93;
t78 = t102 * t191 + t185 * t319;
t281 = t78 * t94;
t79 = t102 * t193 + t187 * t319;
t280 = t79 * t95;
t278 = Ifges(3,4) * t160;
t120 = pkin(1) * mrSges(3,2) + t170 * t175;
t276 = t120 * t28;
t275 = t120 * t29;
t274 = t120 * t30;
t273 = t183 * t28;
t272 = t185 * t29;
t271 = t187 * t30;
t270 = t28 / 0.2e1;
t269 = t29 / 0.2e1;
t268 = t30 / 0.2e1;
t267 = t109 * t151;
t266 = t109 * t154;
t265 = t110 * t152;
t264 = t110 * t155;
t263 = t111 * t153;
t262 = t111 * t156;
t251 = t40 * t304;
t250 = t41 * t303;
t249 = t42 * t302;
t243 = t28 * t278;
t242 = t29 * t278;
t241 = t30 * t278;
t227 = -t142 * t170 + Ifges(2,6);
t67 = t312 * t183 - (t139 * t171 - t227) * t189;
t240 = t164 * t67 * t93;
t226 = -t143 * t170 + Ifges(2,6);
t68 = t313 * t185 - (t140 * t171 - t226) * t191;
t239 = t165 * t68 * t94;
t225 = -t144 * t170 + Ifges(2,6);
t69 = t314 * t187 - (t141 * t171 - t225) * t193;
t238 = t166 * t69 * t95;
t167 = t189 ^ 2;
t86 = t120 * t171 - 0.2e1 * t278 + t315;
t237 = t28 * t86 * t167;
t168 = t191 ^ 2;
t236 = t29 * t86 * t168;
t169 = t193 ^ 2;
t235 = t30 * t86 * t169;
t112 = t183 * pkin(1) + pkin(2) * sin(t161);
t234 = t109 * t112 * t70;
t233 = t109 * t190 * t67;
t113 = t185 * pkin(1) + pkin(2) * sin(t162);
t232 = t110 * t113 * t71;
t231 = t110 * t192 * t68;
t114 = t187 * pkin(1) + pkin(2) * sin(t163);
t230 = t111 * t114 * t72;
t229 = t111 * t194 * t69;
t228 = pkin(2) * t253;
t106 = g(1) * t154 - g(2) * t151;
t218 = g(3) * t190 + t106 * t184;
t107 = g(1) * t155 - g(2) * t152;
t217 = g(3) * t192 + t107 * t186;
t108 = g(1) * t156 - g(2) * t153;
t216 = g(3) * t194 + t108 * t188;
t100 = mrSges(2,1) - t102;
t118 = mrSges(2,2) + t319;
t212 = -t100 * t189 + t118 * t183 - mrSges(1,1);
t211 = -t100 * t191 + t118 * t185 - mrSges(1,1);
t210 = -t100 * t193 + t118 * t187 - mrSges(1,1);
t206 = 0.4e1 * (0.2e1 * t277 + t301) * t171 + 0.2e1 * t101 + 0.4e1 * t255;
t205 = -0.2e1 * t171 * t277 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t255;
t202 = 0.2e1 * qJ(2,1);
t201 = 0.2e1 * qJ(2,2);
t200 = 0.2e1 * qJ(2,3);
t182 = xDDP(1);
t181 = xDDP(2);
t180 = xDDP(3);
t128 = t150 + t323;
t127 = t149 + t323;
t126 = t148 + t323;
t105 = g(1) * t153 + g(2) * t156;
t104 = g(1) * t152 + g(2) * t155;
t103 = g(1) * t151 + g(2) * t154;
t99 = t297 + Ifges(2,3) + Ifges(3,3) + (-t279 + t147) * t253;
t48 = (m(3) * t89 + t194 * t79) * t166;
t47 = (m(3) * t88 + t192 * t78) * t165;
t46 = (m(3) * t87 + t190 * t77) * t164;
t45 = t193 * t86 * t302 + t169 * t76 + (t309 + t300) * qJ(3,1) + t205;
t44 = t191 * t86 * t303 + t168 * t76 + (t309 + t299) * qJ(3,2) + t205;
t43 = t189 * t86 * t304 + t167 * t76 + (t309 + t298) * qJ(3,3) + t205;
t39 = (t194 * t45 + t79 * t89) * t166;
t38 = (t192 * t44 + t78 * t88) * t165;
t37 = (t190 * t43 + t77 * t87) * t164;
t36 = -t66 * t238 + t99 * t263;
t35 = -t65 * t239 + t99 * t265;
t34 = -t64 * t240 + t99 * t267;
t33 = -t63 * t238 + t99 * t262;
t32 = -t62 * t239 + t99 * t264;
t31 = -t61 * t240 + t99 * t266;
t27 = (m(3) * t59 - t66 * t280) * t166;
t26 = (m(3) * t57 - t65 * t281) * t165;
t25 = (m(3) * t55 - t64 * t282) * t164;
t24 = (m(3) * t60 - t63 * t280) * t166;
t23 = (m(3) * t58 - t62 * t281) * t165;
t22 = (m(3) * t56 - t61 * t282) * t164;
t21 = t69 * t263 + (-t45 * t283 + t59 * t79) * t166;
t20 = t68 * t265 + (-t44 * t284 + t57 * t78) * t165;
t19 = t67 * t267 + (-t43 * t285 + t55 * t77) * t164;
t18 = t69 * t262 + (-t45 * t286 + t60 * t79) * t166;
t17 = t68 * t264 + (-t44 * t287 + t58 * t78) * t165;
t16 = t67 * t266 + (-t43 * t288 + t56 * t77) * t164;
t15 = (-0.1e1 / (t159 + (t171 * t193 - t256) * pkin(2)) * t320 + (t213 * t30 + t326) * t30) * t166;
t14 = (-0.1e1 / (t158 + (t171 * t191 - t257) * pkin(2)) * t321 + (t214 * t29 + t327) * t29) * t165;
t13 = (-0.1e1 / (t157 + (t171 * t189 - t258) * pkin(2)) * t322 + (t215 * t28 + t328) * t28) * t164;
t12 = (-t320 + ((-t203 * cos(0.2e1 * t163) - cos(t202) * t204 + (-cos(pkin(5) + t202) - t171) * t228 + t254) * t268 + t316 * t326 + 0.2e1 * (t114 * t75 - t174 * t268) * t174) * t30) * t166;
t11 = (-t321 + ((-t203 * cos(0.2e1 * t162) - cos(t201) * t204 + (-cos(pkin(5) + t201) - t171) * t228 + t254) * t269 + t317 * t327 + 0.2e1 * (t113 * t74 - t173 * t269) * t173) * t29) * t165;
t10 = (-t322 + ((-t203 * cos(0.2e1 * t161) - cos(t200) * t204 + (-cos(t200 + pkin(5)) - t171) * t228 + t254) * t270 + t318 * t328 + 0.2e1 * (t112 * t73 - t172 * t270) * t172) * t28) * t164;
t9 = -t69 * t15 + t99 * t230 + t30 * (0.2e1 * t235 + (t76 * t271 - t311 * t42) * t193 + 0.2e1 * t241 + (mrSges(3,1) * t249 - t274) * t171 + t125 * t249 - t30 * t315) + (t100 * t216 + t105 * t118) * t187 - t193 * (t100 * t105 - t118 * t216);
t8 = -t68 * t14 + t99 * t232 + t29 * (0.2e1 * t236 + (t76 * t272 - t311 * t41) * t191 + 0.2e1 * t242 + (mrSges(3,1) * t250 - t275) * t171 + t125 * t250 - t29 * t315) + (t100 * t217 + t104 * t118) * t185 - t191 * (t104 * t100 - t118 * t217);
t7 = -t67 * t13 + t99 * t234 + t28 * (0.2e1 * t237 + (t76 * t273 - t311 * t40) * t189 + 0.2e1 * t243 + (mrSges(3,1) * t251 - t276) * t171 + t125 * t251 - t28 * t315) + (t100 * t218 + t103 * t118) * t183 - t189 * (t100 * t103 - t118 * t218);
t6 = -t30 ^ 2 * t150 - t79 * t15 + (-g(3) * t188 + t108 * t194 - t12) * m(3) + (-t102 * t187 + t193 * t319) * t75 * t306;
t5 = -t29 ^ 2 * t149 - t78 * t14 + (-g(3) * t186 + t107 * t192 - t11) * m(3) + (-t102 * t185 + t191 * t319) * t74 * t307;
t4 = -t28 ^ 2 * t148 - t77 * t13 + (-g(3) * t184 + t106 * t190 - t10) * m(3) + (-t102 * t183 + t189 * t319) * t73 * t308;
t3 = -t45 * t15 + t69 * t230 - t79 * t12 - t72 * t225 * t187 + t42 * t150 * t306 + (-g(3) * t128 + t108 * t210) * t194 - (g(3) * t210 + t108 * t128) * t188 + (-t206 * t271 * t193 + t315 * t306 + t274 * t305 - 0.4e1 * t235 - 0.4e1 * t241 + (t141 * t187 * t325 + t314 * t193) * t75) * t75;
t2 = -t44 * t14 + t68 * t232 - t78 * t11 - t71 * t226 * t185 + t41 * t149 * t307 + (-g(3) * t127 + t107 * t211) * t192 - (g(3) * t211 + t107 * t127) * t186 + (-t206 * t272 * t191 + t315 * t307 + t275 * t305 - 0.4e1 * t236 - 0.4e1 * t242 + (t140 * t185 * t325 + t313 * t191) * t74) * t74;
t1 = -t43 * t13 + t67 * t234 - t77 * t10 - t70 * t227 * t183 + t40 * t148 * t308 + (-g(3) * t126 + t106 * t212) * t190 - (g(3) * t212 + t106 * t126) * t184 + (-t206 * t273 * t189 + t315 * t308 + t276 * t305 - 0.4e1 * t237 - 0.4e1 * t243 + (t139 * t183 * t325 + t312 * t189) * t73) * t73;
t49 = [t7 * t267 + t8 * t265 + t9 * t263 - g(1) * m(4) + (t36 * t262 + t35 * t264 + t34 * t266) * t181 + (t36 * t263 + t35 * t265 + t34 * t267 + m(4)) * t182 + ((-t21 * t283 + t27 * t59) * t182 + (-t21 * t286 + t27 * t60) * t181 + (t194 * t21 + t27 * t89) * t180 - t3 * t283 + t59 * t6) * t166 + ((-t20 * t284 + t26 * t57) * t182 + (-t20 * t287 + t26 * t58) * t181 + (t192 * t20 + t26 * t88) * t180 - t2 * t284 + t57 * t5) * t165 + ((-t19 * t285 + t25 * t55) * t182 + (-t19 * t288 + t25 * t56) * t181 + (t19 * t190 + t25 * t87) * t180 - t1 * t285 + t55 * t4) * t164; t7 * t266 + t8 * t264 + t9 * t262 - g(2) * m(4) + (t33 * t263 + t32 * t265 + t31 * t267) * t182 + (t33 * t262 + t32 * t264 + t31 * t266 + m(4)) * t181 + ((-t18 * t283 + t24 * t59) * t182 + (-t18 * t286 + t24 * t60) * t181 + (t18 * t194 + t24 * t89) * t180 - t3 * t286 + t60 * t6) * t166 + ((-t17 * t284 + t23 * t57) * t182 + (-t17 * t287 + t23 * t58) * t181 + (t17 * t192 + t23 * t88) * t180 - t2 * t287 + t58 * t5) * t165 + ((-t16 * t285 + t22 * t55) * t182 + (-t16 * t288 + t22 * t56) * t181 + (t16 * t190 + t22 * t87) * t180 - t1 * t288 + t56 * t4) * t164; (-g(3) + t180) * m(4) + ((t153 * t229 - t39 * t283 + t48 * t59) * t182 + (t156 * t229 - t39 * t286 + t48 * t60) * t181 + (t194 * t39 + t48 * t89) * t180 + t194 * t3 + t89 * t6) * t166 + ((t152 * t231 - t38 * t284 + t47 * t57) * t182 + (t155 * t231 - t38 * t287 + t47 * t58) * t181 + (t192 * t38 + t47 * t88) * t180 + t192 * t2 + t88 * t5) * t165 + ((t151 * t233 - t37 * t285 + t46 * t55) * t182 + (t154 * t233 - t37 * t288 + t46 * t56) * t181 + (t190 * t37 + t46 * t87) * t180 + t190 * t1 + t87 * t4) * t164;];
tauX  = t49;
