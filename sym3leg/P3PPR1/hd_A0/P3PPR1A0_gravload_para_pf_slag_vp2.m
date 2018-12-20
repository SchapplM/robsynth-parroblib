% Calculate Gravitation load for parallel robot
% P3PPR1A0
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
%   mass of all robot links (including platform)
% mrSges [3x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:28
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3PPR1A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:28:12
% EndTime: 2018-12-20 17:28:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (66->40), mult. (116->68), div. (0->0), fcn. (86->8), ass. (0->34)
t190 = koppelP(1,1);
t189 = koppelP(2,1);
t188 = koppelP(3,1);
t187 = koppelP(1,2);
t186 = koppelP(2,2);
t185 = koppelP(3,2);
t184 = mrSges(3,1);
t183 = mrSges(3,2);
t182 = xP(3);
t181 = m(1) + m(2);
t180 = legFrame(1,3);
t179 = legFrame(2,3);
t178 = legFrame(3,3);
t177 = cos(t182);
t176 = sin(t182);
t175 = cos(t180);
t174 = cos(t179);
t173 = cos(t178);
t172 = sin(t180);
t171 = sin(t179);
t170 = sin(t178);
t169 = t175 * g(1) + t172 * g(2);
t168 = t174 * g(1) + t171 * g(2);
t167 = t173 * g(1) + t170 * g(2);
t166 = -t172 * g(1) + t175 * g(2);
t165 = -t171 * g(1) + t174 * g(2);
t164 = -t170 * g(1) + t173 * g(2);
t163 = -t176 * t187 + t177 * t190;
t162 = -t176 * t186 + t177 * t189;
t161 = -t176 * t185 + t177 * t188;
t160 = -t176 * t190 - t177 * t187;
t159 = -t176 * t189 - t177 * t186;
t158 = -t176 * t188 - t177 * t185;
t1 = [-g(1) * m(3) + (t170 * t164 + t171 * t165 + t172 * t166) * t181 + (-t167 * t173 - t168 * t174 - t169 * t175) * m(2); -g(2) * m(3) + (-t173 * t164 - t174 * t165 - t175 * t166) * t181 + (-t167 * t170 - t168 * t171 - t169 * t172) * m(2); -(-g(1) * t184 - g(2) * t183) * t176 + t177 * (g(1) * t183 - g(2) * t184) + (-(-t160 * t172 + t163 * t175) * t166 - (-t159 * t171 + t162 * t174) * t165 - (-t158 * t170 + t161 * t173) * t164) * t181 + (-(t160 * t175 + t163 * t172) * t169 - (t159 * t174 + t162 * t171) * t168 - (t158 * t173 + t161 * t170) * t167) * m(2);];
taugX  = t1;
