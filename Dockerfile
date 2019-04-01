FROM trzeci/emscripten
COPY . /src 
RUN npm install
RUN npm run build:test
RUN npm test
# RUN npm run build
# RUN npm start