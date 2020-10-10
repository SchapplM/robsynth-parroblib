% Calculate inertia matrix for parallel robot
% P3RPRR1G1A0
% Use Code from Maple symbolic Code Generation
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:14
% EndTime: 2020-03-09 21:23:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (900->99), mult. (866->151), div. (108->4), fcn. (528->29), ass. (0->81)
t212 = sin(pkin(7));
t252 = pkin(1) * t212;
t251 = m(3) * pkin(2) + mrSges(2,1);
t208 = legFrame(3,3) + qJ(1,3);
t205 = pkin(7) + t208;
t202 = qJ(3,3) + t205;
t196 = sin(t202);
t185 = -pkin(1) * sin(t208) - pkin(2) * sin(t205) - pkin(3) * t196;
t217 = sin(qJ(3,3));
t193 = 0.1e1 / (pkin(1) * sin(pkin(7) + qJ(3,3)) + t217 * pkin(2));
t250 = t185 * t193;
t225 = 0.1e1 / pkin(3);
t249 = t185 * t225;
t209 = legFrame(2,3) + qJ(1,2);
t206 = pkin(7) + t209;
t203 = qJ(3,2) + t206;
t197 = sin(t203);
t186 = -pkin(1) * sin(t209) - pkin(2) * sin(t206) - pkin(3) * t197;
t218 = sin(qJ(3,2));
t194 = 0.1e1 / (pkin(1) * sin(pkin(7) + qJ(3,2)) + t218 * pkin(2));
t248 = t186 * t194;
t247 = t186 * t225;
t210 = legFrame(1,3) + qJ(1,1);
t207 = pkin(7) + t210;
t204 = qJ(3,1) + t207;
t198 = sin(t204);
t187 = -pkin(1) * sin(t210) - pkin(2) * sin(t207) - pkin(3) * t198;
t219 = sin(qJ(3,1));
t195 = 0.1e1 / (pkin(1) * sin(pkin(7) + qJ(3,1)) + t219 * pkin(2));
t246 = t187 * t195;
t245 = t187 * t225;
t199 = cos(t202);
t188 = -pkin(1) * cos(t208) - pkin(2) * cos(t205) - pkin(3) * t199;
t244 = t188 * t193;
t243 = t188 * t225;
t200 = cos(t203);
t189 = -pkin(1) * cos(t209) - pkin(2) * cos(t206) - pkin(3) * t200;
t242 = t189 * t194;
t241 = t189 * t225;
t201 = cos(t204);
t190 = -pkin(1) * cos(t210) - pkin(2) * cos(t207) - pkin(3) * t201;
t240 = t190 * t195;
t239 = t190 * t225;
t238 = t193 * t196;
t237 = t193 * t199;
t236 = t194 * t197;
t235 = t194 * t200;
t234 = t195 * t198;
t233 = t195 * t201;
t232 = 0.2e1 * pkin(1);
t231 = 0.2e1 * pkin(2);
t213 = cos(pkin(7));
t230 = pkin(1) * t213 + pkin(2);
t229 = Ifges(1,3) + Ifges(2,3) + Ifges(3,3) + (m(2) + m(3)) * pkin(1) ^ 2 + m(3) * pkin(2) ^ 2;
t220 = cos(qJ(3,3));
t228 = t220 * mrSges(3,1) - mrSges(3,2) * t217;
t221 = cos(qJ(3,2));
t227 = t221 * mrSges(3,1) - mrSges(3,2) * t218;
t222 = cos(qJ(3,1));
t226 = t222 * mrSges(3,1) - mrSges(3,2) * t219;
t192 = t230 * mrSges(3,1) - mrSges(3,2) * t252;
t191 = -mrSges(3,1) * t252 - t230 * mrSges(3,2);
t184 = t191 * t219 + t192 * t222 + Ifges(3,3);
t183 = t191 * t218 + t192 * t221 + Ifges(3,3);
t182 = t191 * t217 + t192 * t220 + Ifges(3,3);
t181 = t226 * t231 + ((t226 + t251) * t213 - (mrSges(3,1) * t219 + t222 * mrSges(3,2) + mrSges(2,2)) * t212) * t232 + t229;
t180 = t227 * t231 + ((t227 + t251) * t213 - (mrSges(3,1) * t218 + t221 * mrSges(3,2) + mrSges(2,2)) * t212) * t232 + t229;
t179 = t228 * t231 + ((t228 + t251) * t213 - (mrSges(3,1) * t217 + t220 * mrSges(3,2) + mrSges(2,2)) * t212) * t232 + t229;
t178 = (Ifges(3,3) * t239 + t184 * t201) * t195;
t177 = (Ifges(3,3) * t241 + t183 * t200) * t194;
t176 = (Ifges(3,3) * t243 + t182 * t199) * t193;
t175 = (Ifges(3,3) * t245 + t184 * t198) * t195;
t174 = (Ifges(3,3) * t247 + t183 * t197) * t194;
t173 = (Ifges(3,3) * t249 + t182 * t196) * t193;
t172 = (t181 * t201 + t184 * t239) * t195;
t171 = (t180 * t200 + t183 * t241) * t194;
t170 = (t179 * t199 + t182 * t243) * t193;
t169 = (t181 * t198 + t184 * t245) * t195;
t168 = (t180 * t197 + t183 * t247) * t194;
t167 = (t179 * t196 + t182 * t249) * t193;
t1 = [t170 * t237 + t171 * t235 + t172 * t233 + m(4) + (t176 * t244 + t177 * t242 + t178 * t240) * t225, t170 * t238 + t171 * t236 + t172 * t234 + (t176 * t250 + t177 * t248 + t178 * t246) * t225, 0; t167 * t237 + t168 * t235 + t169 * t233 + (t173 * t244 + t174 * t242 + t175 * t240) * t225, t167 * t238 + t168 * t236 + t169 * t234 + m(4) + (t173 * t250 + t174 * t248 + t175 * t246) * t225, 0; 0, 0, 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
