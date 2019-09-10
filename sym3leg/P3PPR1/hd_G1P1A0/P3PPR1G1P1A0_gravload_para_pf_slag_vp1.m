% Calculate Gravitation load for parallel robot
% P3PPR1G1P1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PPR1G1P1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:30
% EndTime: 2019-05-03 14:37:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (66->40), mult. (117->69), div. (0->0), fcn. (86->8), ass. (0->34)
t191 = koppelP(1,1);
t190 = koppelP(2,1);
t189 = koppelP(3,1);
t188 = koppelP(1,2);
t187 = koppelP(2,2);
t186 = koppelP(3,2);
t185 = rSges(3,1);
t184 = rSges(3,2);
t183 = xP(3);
t182 = m(1) + m(2);
t181 = legFrame(1,3);
t180 = legFrame(2,3);
t179 = legFrame(3,3);
t178 = cos(t183);
t177 = sin(t183);
t176 = cos(t181);
t175 = cos(t180);
t174 = cos(t179);
t173 = sin(t181);
t172 = sin(t180);
t171 = sin(t179);
t170 = t176 * g(1) + t173 * g(2);
t169 = t175 * g(1) + t172 * g(2);
t168 = t174 * g(1) + t171 * g(2);
t167 = -t173 * g(1) + t176 * g(2);
t166 = -t172 * g(1) + t175 * g(2);
t165 = -t171 * g(1) + t174 * g(2);
t164 = -t177 * t188 + t178 * t191;
t163 = -t177 * t187 + t178 * t190;
t162 = -t177 * t186 + t178 * t189;
t161 = -t177 * t191 - t178 * t188;
t160 = -t177 * t190 - t178 * t187;
t159 = -t177 * t189 - t178 * t186;
t1 = [-m(3) * g(1) + (t171 * t165 + t172 * t166 + t173 * t167) * t182 + (-t168 * t174 - t169 * t175 - t170 * t176) * m(2); -m(3) * g(2) + (-t174 * t165 - t175 * t166 - t176 * t167) * t182 + (-t168 * t171 - t169 * t172 - t170 * t173) * m(2); m(3) * ((g(1) * t185 + g(2) * t184) * t177 + (g(1) * t184 - g(2) * t185) * t178) + (-(-t161 * t173 + t164 * t176) * t167 - (-t160 * t172 + t163 * t175) * t166 - (-t159 * t171 + t162 * t174) * t165) * t182 + (-(t161 * t176 + t164 * t173) * t170 - (t160 * t175 + t163 * t172) * t169 - (t159 * t174 + t162 * t171) * t168) * m(2);];
taugX  = t1;
