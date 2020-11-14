% Return the minimum parameter vector for
% P3RRPRR8V1G2A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [13x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P3RRPRR8V1G2A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(3,3)}
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_minimal_parameter_para: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V1G2A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_minimal_parameter_para: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
t1238 = Ifges(3,1) - Ifges(3,2);
t1232 = sin(pkin(5));
t1233 = cos(pkin(5));
t1237 = t1232 * t1233;
t1230 = t1232 ^ 2;
t1231 = t1233 ^ 2;
t1236 = t1231 - t1230;
t1235 = Ifges(3,4) * t1237;
t1234 = t1233 * mrSges(3,1) - t1232 * mrSges(3,2);
t1 = [t1230 * Ifges(3,1) + t1231 * Ifges(3,2) + Ifges(2,2) + Ifges(1,3) + 0.2e1 * t1235; mrSges(1,1); mrSges(1,2) - mrSges(2,3); t1238 * t1236 + Ifges(2,1) - Ifges(2,2) - 0.4e1 * t1235; t1236 * Ifges(3,4) + t1238 * t1237 + Ifges(2,4); t1233 * Ifges(3,5) - t1232 * Ifges(3,6) + Ifges(2,5); t1232 * Ifges(3,5) + t1233 * Ifges(3,6) + Ifges(2,6); 0.2e1 * pkin(1) * t1234 + Ifges(2,3) + Ifges(3,3); mrSges(2,1) + t1234; t1232 * mrSges(3,1) + t1233 * mrSges(3,2) + mrSges(2,2); mrSges(3,3); m(3); m(4);];
MPV  = t1;
