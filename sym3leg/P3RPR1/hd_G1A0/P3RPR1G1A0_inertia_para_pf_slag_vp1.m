% Calculate inertia matrix for parallel robot
% P3RPR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPR1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:04
% EndTime: 2019-05-03 14:58:04
% DurationCPUTime: 0.34s
% Computational Cost: add. (1219->121), mult. (1743->207), div. (216->3), fcn. (1304->14), ass. (0->92)
t257 = 2 * rSges(2,3);
t234 = pkin(1) + rSges(2,1);
t256 = m(2) * t234;
t239 = 0.1e1 / qJ(2,3);
t255 = m(2) * t239;
t240 = 0.1e1 / qJ(2,2);
t254 = m(2) * t240;
t241 = 0.1e1 / qJ(2,1);
t253 = m(2) * t241;
t252 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,2) + Icges(1,3);
t251 = rSges(2,3) ^ 2 + pkin(1) ^ 2 + (2 * pkin(1) + rSges(2,1)) * rSges(2,1);
t247 = koppelP(1,1);
t246 = koppelP(2,1);
t245 = koppelP(3,1);
t244 = koppelP(1,2);
t243 = koppelP(2,2);
t242 = koppelP(3,2);
t238 = rSges(3,1);
t237 = rSges(3,2);
t236 = xP(3);
t235 = pkin(1) + pkin(2);
t233 = cos(qJ(1,1));
t232 = cos(qJ(1,2));
t231 = cos(qJ(1,3));
t230 = sin(qJ(1,1));
t229 = sin(qJ(1,2));
t228 = sin(qJ(1,3));
t227 = legFrame(1,3);
t226 = legFrame(2,3);
t225 = legFrame(3,3);
t223 = cos(t236);
t222 = sin(t236);
t221 = cos(t227);
t220 = cos(t226);
t219 = cos(t225);
t218 = sin(t227);
t217 = sin(t226);
t216 = sin(t225);
t214 = t230 * qJ(2,1) + t235 * t233;
t213 = t229 * qJ(2,2) + t235 * t232;
t212 = t228 * qJ(2,3) + t235 * t231;
t211 = -t233 * qJ(2,1) + t230 * t235;
t210 = -t232 * qJ(2,2) + t229 * t235;
t209 = -t231 * qJ(2,3) + t228 * t235;
t208 = -t222 * t244 + t223 * t247;
t207 = -t222 * t243 + t223 * t246;
t206 = -t222 * t242 + t223 * t245;
t205 = -t222 * t247 - t223 * t244;
t204 = -t222 * t246 - t223 * t243;
t203 = -t222 * t245 - t223 * t242;
t202 = m(3) * (-t222 * t237 + t223 * t238);
t201 = m(3) * (-t222 * t238 - t223 * t237);
t200 = -t218 * t230 + t221 * t233;
t199 = t218 * t233 + t221 * t230;
t198 = -t217 * t229 + t220 * t232;
t197 = t217 * t232 + t220 * t229;
t196 = -t216 * t228 + t219 * t231;
t195 = t216 * t231 + t219 * t228;
t194 = ((t257 + qJ(2,1)) * qJ(2,1) + t251) * m(2) + t252;
t193 = ((t257 + qJ(2,2)) * qJ(2,2) + t251) * m(2) + t252;
t192 = ((t257 + qJ(2,3)) * qJ(2,3) + t251) * m(2) + t252;
t191 = -t218 * t211 + t214 * t221;
t190 = -t217 * t210 + t213 * t220;
t189 = -t216 * t209 + t212 * t219;
t188 = t211 * t221 + t218 * t214;
t187 = t210 * t220 + t217 * t213;
t186 = t209 * t219 + t216 * t212;
t185 = (-t200 * t234 + t191) * t253;
t184 = (-t199 * t234 + t188) * t253;
t183 = (-t198 * t234 + t190) * t254;
t182 = (-t197 * t234 + t187) * t254;
t181 = (-t196 * t234 + t189) * t255;
t180 = (-t195 * t234 + t186) * t255;
t179 = (t199 * t208 + t200 * t205) * t241;
t178 = (t197 * t207 + t198 * t204) * t240;
t177 = (t195 * t206 + t196 * t203) * t239;
t176 = (-t191 * t256 + t194 * t200) * t241;
t175 = (-t188 * t256 + t194 * t199) * t241;
t174 = (-t190 * t256 + t193 * t198) * t240;
t173 = (-t187 * t256 + t193 * t197) * t240;
t172 = (-t189 * t256 + t192 * t196) * t239;
t171 = (-t186 * t256 + t192 * t195) * t239;
t170 = (t188 * t208 + t191 * t205) * t241;
t169 = (t187 * t207 + t190 * t204) * t240;
t168 = (t186 * t206 + t189 * t203) * t239;
t167 = (-t179 * t234 + t170) * m(2);
t166 = (-t178 * t234 + t169) * m(2);
t165 = (-t177 * t234 + t168) * m(2);
t164 = -t170 * t256 + t179 * t194;
t163 = -t169 * t256 + t178 * t193;
t162 = -t168 * t256 + t177 * t192;
t1 = [m(3) + (t176 * t200 + t185 * t191) * t241 + (t174 * t198 + t183 * t190) * t240 + (t172 * t196 + t181 * t189) * t239, (t176 * t199 + t185 * t188) * t241 + (t174 * t197 + t183 * t187) * t240 + (t172 * t195 + t181 * t186) * t239, t181 * t168 + t183 * t169 + t185 * t170 + t172 * t177 + t174 * t178 + t176 * t179 + t201; (t175 * t200 + t184 * t191) * t241 + (t173 * t198 + t182 * t190) * t240 + (t171 * t196 + t180 * t189) * t239, m(3) + (t175 * t199 + t184 * t188) * t241 + (t173 * t197 + t182 * t187) * t240 + (t171 * t195 + t180 * t186) * t239, t180 * t168 + t182 * t169 + t184 * t170 + t171 * t177 + t173 * t178 + t175 * t179 + t202; t201 + (t164 * t200 + t167 * t191) * t241 + (t163 * t198 + t166 * t190) * t240 + (t162 * t196 + t165 * t189) * t239, t202 + (t164 * t199 + t167 * t188) * t241 + (t163 * t197 + t166 * t187) * t240 + (t162 * t195 + t165 * t186) * t239, t164 * t179 + t167 * t170 + t163 * t178 + t166 * t169 + t162 * t177 + t165 * t168 + Icges(3,3) + m(3) * (t237 ^ 2 + t238 ^ 2);];
MX  = t1;
